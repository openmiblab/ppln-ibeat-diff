"""
Automatic download of DIXON data from XNAT.
"""
import logging

from miblab import pipe
from miblab_data.xnat import download_series
from tqdm import tqdm


PIPELINE = 'dixon'

SIEMENS = {
    "series_description": [
        "DTI_kidneys_cor-oblique_fb_3mm_23_SLICES_DFC", 
        "DTI_kidneys_cor-oblique_fb_3mm_23_SLICES_DFC_MIX",
        "DTI_kidneys_cor-oblique_fb_DFC",
        "DTI_kidneys_cor-oblique_fb_DFC_MIX",
    ]
} 
TURKU_PHILIPS = {
    "series_description": [
        'DTI-kidney-coronal-oblique-FB',
    ]
}  
TURKU_GE = {
    "series_description": [
        "DTI_kidneys_cor-oblique_fb", 
    ]
} 
TURKU_GE_SETUP = {
    "series_description": [

    ]
} 
SHEFFIELD_SETUP = {
    "series_description": [
        # 'T2star_map_kidneys_cor-oblique_mbh'
    ]
}
SHEFFIELD_PATIENTS = {
    "series_description": [
        # Philips data
        'DTI TEST 2 NSA2',
        'Reg - DTI_3b_RT 20 dir',
        # GE data
        'DTI_kidneys_cor-oblique_fb',
    ]
} 
BARI = {
    "series_description": [
        'DTI TEST 4 NSA1 b100 add',
        'Reg - DTI_3b_RT 20 dir',
    ]
} 
LEEDS = {
    # "parameters/sequence": ["*fl3d2"]
    'DTI_kidneys_cor-oblique_fb',
}
LEEDS_SETUP = {
    "parameters/sequence": ["*fl3d2"],
    "frames": [144]
} 



DOWNLOAD = {
    'leeds_patients':{
        'project_id': "BEAt-DKD-WP4-Leeds",
        'subject_label':"Leeds_Patients",
        'attr': LEEDS      
    },
    'leeds_volunteers':{
        'project_id': "BEAt-DKD-WP4-Leeds",
        'subject_label':"Leeds_volunteer_repeatability_study",
        'attr': LEEDS      
    },
    'leeds_setup':{
        'project_id': "BEAt-DKD-WP4-Leeds",
        'subject_label':"Leeds_setup_scans",
        'attr': LEEDS_SETUP      
    },
    'bari_patients':{
        'project_id': "BEAt-DKD-WP4-Bari",
        'subject_label':"Bari_Patients",
        'attr': BARI      
    },
    'bari_volunteers':{
        'project_id': "BEAt-DKD-WP4-Bari",
        'subject_label':"Bari_Volunteers_Repeatability",
        'attr': BARI      
    },
    'sheffield_patients':{
        'project_id': "BEAt-DKD-WP4-Sheffield",
        'attr': SHEFFIELD_PATIENTS      
    },
    'sheffield_setup':{
        'project_id': "ibeat_setup",
        'attr': SHEFFIELD_SETUP      
    },
    'turku_ge_patients':{
        'project_id': "BEAt-DKD-WP4-Turku",
        'subject_label':"Turku_Patients_GE",
        'attr': TURKU_GE       
    },
    'turku_ge_repeatability':{
        'project_id': "BEAt-DKD-WP4-Turku",
        'subject_label': "Turku_Volunteers_GE_Repeatability",
        'attr': TURKU_GE     
    },
    'turku_ge_setup':{
        'project_id': "BEAt-DKD-WP4-Turku",
        'subject_label': "Turku_GE_Setup_Tests",
        'attr': TURKU_GE_SETUP     
    },
    'turku_philips_patients':{
        'project_id': "BEAt-DKD-WP4-Turku",
        'subject_label': "Turku_Patients_Philips",
        'attr': TURKU_PHILIPS      
    },
    'turku_philips_repeatability':{
        'project_id': "BEAt-DKD-WP4-Turku",
        'subject_label': "Turku_volunteer_repeatability_study",
        'attr': TURKU_PHILIPS    
    },
    'bordeaux_patients_baseline':{
        'project_id': "BEAt-DKD-WP4-Bordeaux",
        'subject_label': "Bordeaux_Patients_Baseline",
        'attr': SIEMENS      
    },
    'bordeaux_volunteers':{
        'project_id': "BEAt-DKD-WP4-Bordeaux",
        'subject_label': "Bordeaux_Volunteers_Repeatability_Baseline",
        'attr': SIEMENS      
    },
    'bordeaux_patients_followup':{
        'project_id': "BEAt-DKD-WP4-Bordeaux",
        'subject_label': "Bordeaux_Patients_Followup",
        'attr': SIEMENS     
    },
    'exeter_patients_baseline':{
        'project_id': "BEAt-DKD-WP4-Exeter",
        'subject_label': "Exeter_Patients_Baseline",
        'attr': SIEMENS      
    },
    'exeter_patients_followup':{
        'project_id': "BEAt-DKD-WP4-Exeter",
        'subject_label': "Exeter_Patients_Followup",
        'attr': SIEMENS      
    },
    'exeter_volunteers':{
        'project_id': "BEAt-DKD-WP4-Exeter",
        'subject_label': "Exeter_Volunteer",
        'attr': SIEMENS     
    },
    'exeter_setup':{
        'project_id': "BEAt-DKD-WP4-Exeter",
        'subject_label': "Exeter_setup_scans",
        'attr': SIEMENS       
    },
}



def run(build, dir_output):

    logging.info("Stage 1 --- Downloading data ---")
 
    n_max=None

    for group, props in tqdm(DOWNLOAD.items(), desc='Downloading..'):
        try:
            download_series(
                xnat_url="https://qib.shef.ac.uk",
                output_dir=dir_output,
                log=True,
                n_max=n_max,
                **props
            )
            logging.info(f"Finished downloading {group}.")
        except:
            logging.exception(f"Error downloading {group}.")



if __name__ == '__main__':

    build = r"C:\Users\eic20eh\Documents\Data\iBEAt_Build"
    pipe.run_stage(run, build, PIPELINE, __file__)