import pickle
import logging

def save_obj(obj, name):
    '''
    @param obj: Dictionary of k-mers
    @param name: Name of the output file
    @return: Create a pickle file with the content the obj
    '''
    with open(name + '.pkl', 'wb') as f:
        logging.info('Saving pickle file with name: ' + name + '.pkl')
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    '''
    @param name: Name of the pickle file to read
    @return: Load in a dictionary the content from the pickle file
    '''
    with open(name + '.pkl', 'rb') as f:
        logging.info('Loading pickle file with name: ' + name + '.pkl')
        return pickle.load(f)
