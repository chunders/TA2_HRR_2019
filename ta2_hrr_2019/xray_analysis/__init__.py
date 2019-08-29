try:
    from . import paths
except ImportError:
    print('You have not specified the path to the data.')
    print('Create a paths.py file by copying paths-template.py and editing it.')
    raise
