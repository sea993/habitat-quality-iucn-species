from pathlib import Path
import argparse
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import r2_score
from sklearn.neighbors import KDTree
from sklearn.cluster import MiniBatchKMeans
import shap
import joblib
import warnings
import time
import os
import gc
import shutil
from tqdm import tqdm

warnings.filterwarnings('ignore')


# ============================================
# 1. 基础配置与数据加载
# ============================================

def load_all_data_once(file_path):
    print(f"📂 加载数据文件: {file_path}")
    start_time = time.time()
    try:
        data = pd.read_parquet(file_path) if Path(file_path).suffix == '.parquet' else pd.read_excel(file_path, engine='openpyxl')
        if 'X' in data.columns and 'Y' in data.columns:
            coords = data[['X', 'Y']]
        elif 'longitude' in data.columns and 'latitude' in data.columns:
            coords = data[['longitude', 'latitude']]
            coords.columns = ['X', 'Y']
        else:
            print("⚠ 警告: 未找到坐标列，使用索引生成坐标")
            coords = pd.DataFrame({'X': np.arange(len(data)), 'Y': np.zeros(len(data))})
    except Exception as e:
        print(f"❌ 数据加载失败: {e}")
        raise
    end_time = time.time()
    print(f"✓ 数据加载完成! 耗时: {(end_time - start_time):.2f} 秒, 形状: {data.shape}")
    data.columns = [c.strip() for c in data.columns]
    return data, coords


def calculate_global_stats(data, years):
    print("\n📊 计算全局归一化参数 (2000-2020)...")
    dynamic_features = ['NDVI', 'GDP', 'PRE', 'POP', 'CE', 'HFI', 'TMP', 'CSDI', 'GSL', 'WSDI']
    static_features = ['DEM', 'SLOPE']
    target_base = 'HQ'
    stats = {}

    for feat in dynamic_features:
        all_values = []
        for year in years:
            col_name = f"{feat}{year}"
            if col_name in data.columns: all_values.append(data[col_name])
        if all_values:
            combined = pd.concat(all_values)
            stats[feat] = {'min': combined.min(), 'max': combined.max()}

    for feat in static_features:
        if feat in data.columns:
            stats[feat] = {'min': data[feat].min(), 'max': data[feat].max()}

    hq_values = []
    for year in years:
        col_name = f"{target_base}{year}"
        if col_name in data.columns: hq_values.append(data[col_name])
    if hq_values:
        combined_hq = pd.concat(hq_values)
        stats['HQ'] = {'min': combined_hq.min(), 'max': combined_hq.max()}
        print(f"  Target (HQ) 全局范围: [{stats['HQ']['min']:.4f}, {stats['HQ']['max']:.4f}]")
    return stats


def normalize_data(df, stats, year, is_target=False, target_name='HQ'):
    df_norm = df.copy()
    if is_target:
        min_val = stats[target_name]['min']
        max_val = stats[target_name]['max']
        if max_val - min_val == 0: return df_norm
        return (df_norm - min_val) / (max_val - min_val)
    else:
        for col in df.columns:
            base_name = None
            if col in stats:
                base_name = col
            else:
                for key in stats.keys():
                    if col.startswith(key) and col != key:
                        if col.replace(year, "") == key:
                            base_name = key
                            break

            if base_name and base_name in stats:
                min_val = stats[base_name]['min']
                max_val = stats[base_name]['max']
                if max_val - min_val != 0:
                    df_norm[col] = (df[col] - min_val) / (max_val - min_val)
                else:
                    df_norm[col] = 0
        return df_norm


def inverse_normalize_target(y_norm, stats, target_name='HQ'):
    min_val = stats[target_name]['min']
    max_val = stats[target_name]['max']
    return y_norm * (max_val - min_val) + min_val


def extract_year_data(data, coords, year):
    target_col = f'HQ{year}'
    potential_features = [
        f'NDVI{year}', f'GDP{year}', f'PRE{year}', f'POP{year}',
        f'CE{year}', f'HFI{year}', f'TMP{year}', f'CSDI{year}',
        f'GSL{year}', f'WSDI{year}', 'DEM', 'SLOPE'
    ]
    available_features = [col for col in potential_features if col in data.columns]
    if target_col not in data.columns: raise ValueError(f"找不到目标变量 {target_col}")

    features = data[available_features]
    target = data[target_col]
    valid_idx = ~(features.isna().any(axis=1) | target.isna() | coords.isna().any(axis=1))

    return (features[valid_idx].reset_index(drop=True),
            target[valid_idx].reset_index(drop=True),
            coords[valid_idx].reset_index(drop=True),
            available_features)


# ============================================
# 2. GWRF-SHAP 核心模型 (修改版：支持保存元数据)
# ============================================

class GWRF_SHAP_Synced:
    def __init__(self, bandwidth, rf_params, max_local_models=1000, n_jobs=16, temp_dir=None):
        self.bandwidth = bandwidth
        self.rf_params = rf_params
        self.max_local_models = max_local_models
        self.n_jobs = n_jobs
        self.temp_dir = temp_dir
        self.model_centers = None
        self.model_tree = None
        self.feature_names = None
        self.valid_models = set()
        self.model_indices = None
        self.data_cache_path = os.path.join(temp_dir, 'training_data_cache.joblib')

    @staticmethod
    def _fit_single_model_static(idx, data_cache_path, bandwidth, rf_params, temp_dir):
        try:
            X_np, y_np, coords_np = joblib.load(data_cache_path, mmap_mode='r')
            target_point = coords_np[idx]

            distances = np.sqrt(np.sum((coords_np - target_point) ** 2, axis=1))
            weights = np.exp(-0.5 * (distances / bandwidth) ** 2)

            mask = weights > 1e-4

            if np.sum(mask) < 20:
                top_k_idx = np.argsort(weights)[-50:]
                mask = np.zeros_like(mask, dtype=bool)
                mask[top_k_idx] = True

            X_local = X_np[mask]
            y_local = y_np[mask]
            w_local = weights[mask]

            if len(np.unique(y_local)) < 2:
                return idx, False

            local_params = rf_params.copy()
            local_params['n_estimators'] = min(local_params.get('n_estimators', 100), 60)
            local_params['n_jobs'] = 1

            model = RandomForestRegressor(**local_params)
            model.fit(X_local, y_local, sample_weight=w_local)

            # 保存模型文件
            model_path = os.path.join(temp_dir, f'model_{idx}.pkl')
            joblib.dump(model, model_path, compress=3)

            return idx, True
        except Exception as e:
            # print(f"  ⚠ 模型 {idx} 训练失败: {e}")
            return idx, False

    def fit(self, X, y, coords, feature_names):
        self.feature_names = feature_names

        X_np = X.values.astype(np.float32)
        y_np = y.values.astype(np.float32)
        coords_np = coords.values.astype(np.float32)
        joblib.dump((X_np, y_np, coords_np), self.data_cache_path)

        n_samples = len(coords_np)

        if n_samples > self.max_local_models:
            print(f"  执行空间聚类 ({n_samples} -> {self.max_local_models} 模型中心)...")
            kmeans = MiniBatchKMeans(n_clusters=self.max_local_models,
                                     batch_size=2048,
                                     random_state=42,
                                     n_init='auto').fit(coords_np)
            self.model_centers = kmeans.cluster_centers_
            tree = KDTree(coords_np)
            _, idxs = tree.query(self.model_centers, k=1)
            model_indices = idxs.flatten()
        else:
            self.model_centers = coords_np
            model_indices = np.arange(n_samples)

        self.model_tree = KDTree(self.model_centers)
        self.model_indices = model_indices

        print(f"  开始并行训练 {len(model_indices)} 个局部模型...")
        results = joblib.Parallel(n_jobs=self.n_jobs)(
            joblib.delayed(GWRF_SHAP_Synced._fit_single_model_static)(
                idx, self.data_cache_path, self.bandwidth, self.rf_params, self.temp_dir
            )
            for idx in tqdm(model_indices, desc="GWRF训练", unit="模型", leave=False)
        )
        self.valid_models = {res[0] for res in results if res[1]}
        print(f"  ✓ 成功训练 {len(self.valid_models)}/{len(model_indices)} 个模型")

        # ✅ 新增：保存元数据，供后续IS计算脚本使用
        metadata = {
            'model_centers': self.model_centers,
            'model_indices': self.model_indices,  # 聚类中心对应的训练集行索引
            'valid_models': self.valid_models,
            'feature_names': self.feature_names
        }
        meta_path = os.path.join(self.temp_dir, 'model_metadata.joblib')
        joblib.dump(metadata, meta_path)
        print(f"  📝 模型元数据已保存至: {meta_path}")

    def predict_and_calculate_shap(self, X, coords):
        coords_np = coords.values.astype(np.float32)
        X_np = X.values.astype(np.float32)
        n_samples = len(coords)
        n_features = X.shape[1]

        _, nearest_model_indices = self.model_tree.query(coords_np, k=1)
        nearest_model_indices = nearest_model_indices.flatten()

        predictions = np.zeros(n_samples)
        local_shap_values = np.zeros((n_samples, n_features))
        unique_models = np.unique(nearest_model_indices)

        print(f"  正在计算 SHAP 值 (解释变量方向性贡献)...")

        shap_success_count = 0
        shap_fail_count = 0

        for tree_idx in tqdm(unique_models, desc="SHAP计算", leave=False):
            model_idx = self.model_indices[tree_idx]

            if model_idx not in self.valid_models:
                continue

            model_path = os.path.join(self.temp_dir, f'model_{model_idx}.pkl')
            try:
                if not os.path.exists(model_path):
                    shap_fail_count += 1
                    continue

                model = joblib.load(model_path)

                sample_mask = (nearest_model_indices == tree_idx)
                X_subset = X_np[sample_mask]

                if len(X_subset) > 0:
                    predictions[sample_mask] = model.predict(X_subset)
                    explainer = shap.TreeExplainer(model)
                    shap_vals = explainer.shap_values(X_subset, check_additivity=False)

                    if isinstance(shap_vals, list):
                        shap_vals = shap_vals[0]
                    local_shap_values[sample_mask, :] = shap_vals
                    shap_success_count += len(X_subset)

                del model, explainer
                gc.collect()

            except Exception as e:
                shap_fail_count += 1
                continue

        print(f"  ✓ SHAP计算完成：成功 {shap_success_count} 个样本")
        if shap_success_count != len(X_np):
            raise RuntimeError(f'Incomplete SHAP predictions: {shap_success_count}/{len(X_np)}')
        shap_df = pd.DataFrame(local_shap_values, columns=[f'SHAP_{f}' for f in self.feature_names])
        return predictions, shap_df


# ============================================
# 3. 参数优化
# ============================================

def optimize_rf_params(data, coords, year_list, global_stats, n_jobs=8):
    print("正在优化随机森林基础参数...")
    X_list, y_list = [], []

    for year in year_list[:3]:
        try:
            f, t, _, _ = extract_year_data(data, coords, year)
            f_norm = normalize_data(f, global_stats, year, is_target=False)
            t_norm = normalize_data(t, global_stats, year, is_target=True)
            if len(f_norm) > 2000:
                idx = np.random.choice(len(f_norm), 2000, replace=False)
                X_list.append(f_norm.iloc[idx])
                y_list.append(t_norm.iloc[idx])
            else:
                X_list.append(f_norm)
                y_list.append(t_norm)
        except:
            continue

    if not X_list: return {'n_estimators': 100, 'max_depth': 15}
    X_sample = pd.concat(X_list)
    y_sample = pd.concat(y_list)

    param_grid = {
        'n_estimators': [100, 150],
        'max_depth': [10, 20],
        'min_samples_leaf': [2, 5]
    }
    grid = GridSearchCV(RandomForestRegressor(n_jobs=1), param_grid, cv=3, n_jobs=n_jobs)
    grid.fit(X_sample, y_sample)
    print(f"最优参数: {grid.best_params_}")
    return grid.best_params_


# ============================================
# 4. 单年处理流程 (修改：保存路径)
# ============================================

def process_year(year, data, coords, rf_params, global_stats, csv_save_dir, model_save_base, n_jobs):
    print(f"\n{'=' * 30} 处理 {year} 年 {'=' * 30}")

    try:
        X_raw, y_raw, coords_curr, feature_names = extract_year_data(data, coords, year)
    except Exception as e:
        print(f"跳过 {year} 年: {e}")
        return None

    print("  执行 0-1 归一化...")
    X_norm = normalize_data(X_raw, global_stats, year, is_target=False)
    y_norm = normalize_data(y_raw, global_stats, year, is_target=True)

    indices = np.arange(len(X_norm))
    X_train, X_test, y_train, y_test, c_train, c_test, idx_train, idx_test = train_test_split(
        X_norm, y_norm, coords_curr, indices, test_size=0.3, random_state=42
    )

    coord_range = np.sqrt(
        (c_train['X'].max() - c_train['X'].min()) ** 2 + (c_train['Y'].max() - c_train['Y'].min()) ** 2)
    bandwidth = coord_range * 0.15
    print(f"  使用带宽: {bandwidth:.2f}")

    # ✅ 关键修改：指定持久化保存路径
    model_save_dir = os.path.join(model_save_base, year)
    if not os.path.exists(model_save_dir):
        os.makedirs(model_save_dir)
    print(f"  📂 模型将保存至: {model_save_dir}")

    try:
        n_samples = len(X_train)
        max_models = min(2000, int(n_samples / 50))
        if max_models < 100: max_models = 100

        # 初始化模型，传入持久化路径
        model = GWRF_SHAP_Synced(bandwidth, rf_params, max_local_models=max_models, n_jobs=n_jobs,
                                 temp_dir=model_save_dir)
        model.fit(X_train, y_train, c_train, feature_names)

        print("  正在生成全量 SHAP 解释表...")
        pred_all_norm, shap_all_df = model.predict_and_calculate_shap(X_norm, coords_curr)

        pred_all_orig = inverse_normalize_target(pred_all_norm, global_stats)

        output_df = coords_curr.copy()
        output_df['Observed'] = y_raw
        output_df['Predicted'] = pred_all_orig
        output_df['Residual'] = output_df['Observed'] - output_df['Predicted']
        output_df['Set_Type'] = 'Train'
        output_df.loc[idx_test, 'Set_Type'] = 'Test'

        output_df = pd.concat([output_df, shap_all_df], axis=1)

        csv_path = os.path.join(csv_save_dir, f'GWRF_SHAP_Result_{year}.csv')
        output_df.to_csv(csv_path, index=False)
        print(f"  ✓ 详情 (含SHAP) 已保存至: {csv_path}")

        test_mask = output_df['Set_Type'] == 'Test'
        r2 = r2_score(output_df.loc[test_mask, 'Observed'], output_df.loc[test_mask, 'Predicted'])
        print(f"  测试集 R2: {r2:.4f}")

    except Exception as e:
        print(f"  ⚠ 处理年份 {year} 时发生错误: {e}")
        import traceback
        raise

    # Model files are retained for reuse.
    # finally:
    #     try:
    #         shutil.rmtree(temp_dir)
    #     except: pass

    gc.collect()
    return {'Year': year, 'R2_Test': r2}


# ============================================
# 5. 主程序
# ============================================

def main():
    parser = argparse.ArgumentParser(description="Refit the annual GWRF training configuration.")
    china = Path(__file__).resolve().parents[1]
    parser.add_argument('--input', type=Path, default=china / 'model_inputs/huigui_V2.parquet')
    parser.add_argument('--output-dir', type=Path, default=china.parent.parent / 'review-output/refit')
    parser.add_argument('--jobs', type=int, default=max(1, (os.cpu_count() or 1)-4))
    parser.add_argument('--seed', type=int, default=42, help='Seed for a new refit; historical tuning subsamples were unseeded.')
    args = parser.parse_args()
    np.random.seed(args.seed)
    FILE_PATH, CSV_SAVE_DIR, MODEL_SAVE_BASE = args.input, args.output_dir / 'annual', args.output_dir / 'models'
    os.makedirs(CSV_SAVE_DIR, exist_ok=True)
    os.makedirs(MODEL_SAVE_BASE, exist_ok=True)

    n_cores = os.cpu_count()
    n_jobs = args.jobs
    print(f"GWRF 训练模式 (保存PKL) | 系统核心数: {n_cores}, 使用线程数: {n_jobs}")

    data, coords = load_all_data_once(FILE_PATH)
    years = [str(y) for y in range(2000, 2021)]
    global_stats = calculate_global_stats(data, years)

    rf_params = optimize_rf_params(data, coords, years, global_stats, n_jobs=n_jobs)

    for year in years:
        process_year(year, data, coords, rf_params, global_stats, CSV_SAVE_DIR, MODEL_SAVE_BASE, n_jobs)

    print("\n" + "=" * 50)
    print("所有年份训练及保存完成！")


if __name__ == '__main__':
    main()
