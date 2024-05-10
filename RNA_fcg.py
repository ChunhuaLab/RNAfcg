from cif_read import cif_read
from B_factor_nor import b_normal
from scipy.stats import pearsonr
import pandas as pd
import pickle


def RF_Predice():
    B_pre = open("Predict_B_factor.csv", "w")
    all_futere = pd.read_csv("4LNT_RA_all_feature.csv")      # Change to the target PDB file
    test_all_x = all_futere.drop('labels', axis=1)
    test_all_y = all_futere['labels']

    f1 = open('model_RF_final.pkl', "rb")
    model_rf = pickle.load(f1)
    test_predict = model_rf.predict(test_all_x)
    for i in test_predict:
        B_pre.write(str(i))
        B_pre.write("\n")
    pcc = pearsonr(test_all_y, test_predict)
    print("Coefficient of correlation ：", pcc)


if __name__ == '__main__':
    RF_Predice()