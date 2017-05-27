#!/home/quan/anaconda3/bin/python

import pandas as pd
import pickle

class Data_accessor(object):

    def __init__(self, file_name):
        self.file_name = file_name

    def load(self):
        with open(self.file_name, 'rb') as datafile:
            data = pickle.load(datafile)
            return data

        return None

    def dump(self, output_file_name, data, end=None, start=None):
        assert(data is not None)
        assert(output_file_name is not None)

        file_name_tuple = output_file_name.split('.pkl')
        file_name_prefix = file_name_tuple[0]

        if start == None and end != None:
            file_name_prefix += "_0_" + str(end)
            data = data[:end]
        elif start != None and end == None:
            file_name_prefix += "_" + str(start) + "_END"
            data = data[start:]
        elif start != None and end != None:
            file_name_prefix += "_" + str(start) + "_" + str(end)
            data = data[start:end]
        else:
            file_name_prefix += "_0_END"

        output_file_name = file_name_prefix + ".pkl"
        file_handler = open(output_file_name,"wb")
        pickle.dump(data, file_handler)

        return output_file_name

if __name__ == "__main__":

    #file_name = './train/ml_finalproj_train_vF.pkl'
    #da = Data_accessor(file_name)
    #data = da.load()
    #da.dump('./train/ml_finalproj_train_vF_one_third.pkl', data, start=94108)

    train_file_name = './train/ml_finalproj_train_vF_two_thirds.pkl'
    prod_file_name = './train/ml_finalproj_train_vF_one_third.pkl'
    
    da1 = Data_accessor(train_file_name)
    train_data = da1.load()
    
    da2 = Data_accessor(prod_file_name)
    prod_data = da2.load()

    print(train_data.tail())
    print(len(train_data))

    print(prod_data.head())
    print(len(prod_data))




    


