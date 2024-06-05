import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator
from sklearn.model_selection import GridSearchCV


class Stacking(BaseEstimator):

    '''
    This package allows to combine outputs of different Machine Learning models built with different training series on a test series with the same features.
    
    You have to introduce the following hyperparameters:

    model_params_data = {
    'RandomForest': (RandomForestClassifier(random_state=42), 
                     {'n_estimators': [50, 100], 'max_depth': [None, 10]}, 
                     (X1_train, y1_train)),
    'SVM': (SVC(probability=True, random_state=42), 
            {'C': [0.1, 1, 10], 'kernel': ['linear', 'rbf']}, 
            (X2_train, y2_train)),
    'LogisticRegression': (LogisticRegression(random_state=42), 
                           {'C': [0.1, 1, 10]}, 
                           (X3_train, y3_train)),
    'KNN': (KNeighborsClassifier(), 
            {'n_neighbors': [3, 5, 7]}, 
            (X4_train, y4_train))}

    So, in model_params_data you have to introduce a dictionary of so much models as you want with their hyperparameters to be optimized together along their training series.

    logical_rule contains: OR, AND, Majority
    
    '''
    def __init__(self, models_params_data=None, logical_rule='OR'):
        
        if models_params_data is None:
            raise ValueError("models_params_data must be provided")
            
        self.models_params_data = models_params_data
        self.logical_rule = logical_rule
        self.fitted_models = {}
        
    def fit(self, X, y):

        '''
        X and y are test series over which stacking model is built'''
        
        for name, (model, params, (X_train, y_train)) in self.models_params_data.items():
            grid_search = GridSearchCV(model, params, cv=5, scoring='roc_auc')
            grid_search.fit(X_train, y_train)
            self.fitted_models[name] = grid_search.best_estimator_
        return self
    
    def predict(self, X):
        # Obtener predicciones de cada modelo
        predictions = {}
        for name, model in self.fitted_models.items():
            predictions[name] = model.predict(X)  # Obtener probabilidad de la clase positiva

        # Convertir las predicciones en un DataFrame
        df_predictions = pd.DataFrame(predictions)

        # Aplicar la regla lógica
        if df_predictions.shape[1] < 3:
            raise ValueError('Introduce more than one or two individual predictions')

        elif df_predictions.shape[1] >= 3:

            if self.logical_rule == 'OR':
                pred = np.where(np.any(df_predictions, axis=1), 1, 0)
            elif self.logical_rule == 'AND':
                pred = np.where(np.all(df_predictions, axis=1), 1, 0)
            elif self.logical_rule == 'Majority':
                pred = np.where(np.sum(df_predictions, axis=1) > df_predictions.shape[1] / 2, 1, 0)
        else:
            raise ValueError("Invalid individual prediction combination and/or logical rules.")
        
        return pred

    def get_params(self, deep=True):
        return {'models_params_data': self.models_params_data, 'logical_rule': self.logical_rule}
    
    def set_params(self, **parameters):
        for parameter, value in parameters.items():
            setattr(self, parameter, value)
        return self