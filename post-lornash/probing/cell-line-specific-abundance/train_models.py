# !pip install pandas scikit-learn numpy scipy joblib xgboost lightgbm catboost

import pandas as pd
import numpy as np
from sklearn.model_selection import GroupKFold, GridSearchCV
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, AdaBoostRegressor
from sklearn.svm import SVR
from sklearn.tree import ExtraTreeRegressor
from sklearn.neighbors import KNeighborsRegressor
from xgboost import XGBRegressor
from lightgbm import LGBMRegressor
from catboost import CatBoostRegressor
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr
from joblib import Parallel, delayed
import os

def read_data():

    features = pd.read_csv('data/c2l2.lornash.embeds.csv')
    metadata = pd.read_csv('data/c2l2.transcripts.row_data.csv')
    target_data = pd.read_csv('data/c2l2.mpaqt.abun.csv')
    
    return features, metadata, target_data

def preprocess_data(cell_line, features, metadata, target_data):

    data = pd.merge(features, metadata, on='transcript_id')
    
    target = target_data[['transcript_id', cell_line]].copy()
    data = pd.merge(data, target, on='transcript_id')
    data[cell_line] = np.log(data[cell_line])
    
    X = data[[f'dim{i}' for i in range(1, 129)]].values
    y = data[cell_line].values

    groups = data['chr']
    
    return X, y, groups, data

def cross_validate_fold(X_train, X_test, y_train, y_test, model):

    model.fit(X_train, y_train)

    y_pred_train = model.predict(X_train)
    y_pred_test = model.predict(X_test)
    
    mse_train = mean_squared_error(y_train, y_pred_train)
    pearson_train, pvalue_train = pearsonr(y_train, y_pred_train)
    
    mse_test = mean_squared_error(y_test, y_pred_test)
    pearson_test, pvalue_test = pearsonr(y_test, y_pred_test)
    
    return mse_train, pearson_train, pvalue_train, mse_test, pearson_test, pvalue_test

def custom_train_test_split(data, cell_line, model, param_grid=None):

    train_data = data[~data['chr'].isin(['chr10', 'chr8'])]
    val_data = data[data['chr'] == 'chr10']
    test_data = data[data['chr'] == 'chr8']
    
    X_train = train_data[[f'dim{i}' for i in range(1, 129)]].values
    y_train = train_data[cell_line].values
    
    X_val = val_data[[f'dim{i}' for i in range(1, 129)]].values
    y_val = val_data[cell_line].values
    
    X_test = test_data[[f'dim{i}' for i in range(1, 129)]].values
    y_test = test_data[cell_line].values
    
    if param_grid:
        grid_search = GridSearchCV(model, param_grid, cv=3, n_jobs=-1, scoring='neg_mean_squared_error')
        grid_search.fit(X_train, y_train)
        best_model = grid_search.best_estimator_
        best_params = grid_search.best_params_
        print(f"Best parameters for {model.__class__.__name__}: {best_params}")
    else:
        best_model = model
        best_params = {}
    
    X_train_val = np.concatenate([X_train, X_val], axis=0)
    y_train_val = np.concatenate([y_train, y_val], axis=0)

    best_model.fit(X_train_val, y_train_val)

    mse_train, pearson_train, pvalue_train = evaluate_model(best_model, X_train_val, y_train_val)
    mse_val, pearson_val, pvalue_val = evaluate_model(best_model, X_val, y_val)
    mse_test, pearson_test, pvalue_test = evaluate_model(best_model, X_test, y_test)
    
    return mse_train, pearson_train, pvalue_train, mse_val, pearson_val, pvalue_val, mse_test, pearson_test, pvalue_test, best_params

def evaluate_model(model, X, y):
    y_pred = model.predict(X)
    mse = mean_squared_error(y, y_pred)
    pearson_corr, p_value = pearsonr(y, y_pred)
    return mse, pearson_corr, p_value

def train_and_evaluate_models(cell_line, validation_method="custom-split", n_jobs=-1):

    features, metadata, target_data = read_data()
    X, y, groups, data = preprocess_data(cell_line, features, metadata, target_data)
    
    models = {
        'Linear Regression': (LinearRegression(), None),
        'Ridge Regression': (Ridge(), {'alpha': [0.1, 1, 10]}),
        'Lasso Regression': (Lasso(), {'alpha': [0.01, 0.1, 1]}),
        'Random Forest': (RandomForestRegressor(random_state=42), {
            'n_estimators': [50, 100, 200],
            'max_depth': [3, 5, 10]
        }),
        'Gradient Boosting': (GradientBoostingRegressor(random_state=42), {
            'n_estimators': [50, 100],
            'learning_rate': [0.01, 0.1],
            'max_depth': [3, 5]
        }),
        'AdaBoost': (AdaBoostRegressor(random_state=42), {
            'n_estimators': [50, 100]
        }),
        'Support Vector Regressor': (SVR(), {
            'C': [0.1, 1, 10],
            'kernel': ['linear', 'rbf']
        }),
        'XGBoost': (XGBRegressor(random_state=42), {
            'n_estimators': [50, 100],
            'learning_rate': [0.01, 0.1],
            'max_depth': [3, 5]
        }),
        'LightGBM': (LGBMRegressor(random_state=42), {
            'n_estimators': [50, 100],
            'learning_rate': [0.01, 0.1],
            'max_depth': [3, 5]
        }),
        'CatBoost': (CatBoostRegressor(random_state=42, verbose=0), {
            'iterations': [50, 100],
            'learning_rate': [0.01, 0.1],
            'depth': [3, 5]
        })
    }

    results_df = pd.DataFrame(columns=['cell_line', 'model', 'hyperparameters',
                                       'mse_train', 'pearson_train', 'pvalue_train',
                                       'mse_val', 'pearson_val', 'pvalue_val',
                                       'mse_test', 'pearson_test', 'pvalue_test'])
    
    if validation_method == "custom-split":
        print(f"Running custom train-validation-test split for {cell_line}...")

        for model_name, (model, param_grid) in models.items():
            mse_train, pearson_train, pvalue_train, mse_val, pearson_val, pvalue_val, mse_test, pearson_test, pvalue_test, best_params = custom_train_test_split(
                data, cell_line, model, param_grid)

            print(f"Cell Line: {cell_line} | Model: {model_name} - Train MSE: {mse_train}, Validation MSE: {mse_val}, Test MSE: {mse_test}")
            print(f"Cell Line: {cell_line} | Model: {model_name} - Train Pearson: {pearson_train}, Validation: {pearson_val}, Test: {pearson_test}")
            
            temp_df = pd.DataFrame([{
                'cell_line': cell_line,
                'model': model_name,
                'hyperparameters': best_params,
                'mse_train': mse_train,
                'pearson_train': pearson_train,
                'pvalue_train': pvalue_train,
                'mse_val': mse_val,
                'pearson_val': pearson_val,
                'pvalue_val': pvalue_val,
                'mse_test': mse_test,
                'pearson_test': pearson_test,
                'pvalue_test': pvalue_test
            }])

            results_df = pd.concat([results_df, temp_df], ignore_index=True)
    
    results_file = f'results_{cell_line}.csv'
    results_df.to_csv(results_file, index=False)
    print(f"Results saved to {results_file}")

def run_all_cell_lines_parallel(cell_lines, n_jobs=-1):
    Parallel(n_jobs=n_jobs)(delayed(train_and_evaluate_models)(cell_line, n_jobs=n_jobs) for cell_line in cell_lines)

cell_lines = ['A2058', 'A375tr', 'A549', 'BxPC3', 'C42B', 'Colo320', 'H1299', 'H23', 'H358', 'HCC44', 
              'HCT116', 'HEK293', 'HEPG2', 'HS707A', 'K562', 'LNCaP', 'LS174Ttr', 'MCF7', 'MDA_231u',
              'MDA_453', 'MP2', 'PANC1tr', 'PC3', 'SW480tr', 'SW620', 'ZR_75_1']

run_all_cell_lines_parallel(cell_lines, n_jobs=-1)
