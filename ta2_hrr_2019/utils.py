import os

from .paths import DATA_FOLDER

def shots_sorted(files, ext):
    """Sort a list of shot files into numerical order. Selects only a single extension"""
    l = len(ext)
    ext = ext.lower()
    files = [f for f in files if f.lower().endswith(ext)]
    files.sort(key=lambda f: int(f[4:-l]))
    return files

def setup_mirage_analysis():
    import mirage_analysis
    from mirage_analysis.loader import GCamDataLoader, FallibleGCamDataLoader, ImageDataLoader
    mirage_analysis.configure(os.path.join(DATA_FOLDER, 'MIRAGE'))
    mirage_analysis.register_data_loader('ESpec', FallibleGCamDataLoader)
    mirage_analysis.register_data_loader('HighESpec', ImageDataLoader)
    mirage_analysis.register_data_loader('HASOFF', GCamDataLoader)
    mirage_analysis.register_data_loader('HASONF', GCamDataLoader)
    mirage_analysis.register_data_loader('PreCompFF', GCamDataLoader)
    mirage_analysis.register_data_loader('PreCompNF', GCamDataLoader)
    mirage_analysis.register_data_loader('XRay', ImageDataLoader)
    mirage_analysis.register_data_loader('EProfile', GCamDataLoader)
    mirage_analysis.register_data_loader('Probe_Interferometry', ImageDataLoader)
    mirage_analysis.register_data_loader('SideSpec', FallibleGCamDataLoader)
