1. Download and Installation 
pip install .
#Successfully installed sctcm-0.1.0


2. Structure of sctcm
sctcm/                              
├── __init__.py                    
├── sc/                          
│   ├── __init__.py                
│   ├── io.py                     
│   ├── mouse2human.py           
│   ├── qc.py                   
│   ├── dimred.py              
│   ├── cluster.py            
│   ├── marker4anno.py       
│   ├── AUCell4anno.py      
│   ├── recode4anno.py     
│   ├── deg.py            
│   ├── subtypeProc.py            
│   └── marker4subanno.py        
├── tcm/                        
│   ├── __init__.py
│   ├── ingredient2target.py    
│   ├── herb2target.py         
│   ├── formula2target.py     
│   ├── ing2targetall.py     
│   ├── herb2targetall.py   
│   ├── form2targetall.py        
│   ├── herb2nonnested.py       
│   └── form2nonnested.py      
├── ChEMBL/
│   ├── __init__.py
│   ├── parseChEMBLfiltered.py     
│   └── parseChEMBLall.py         
├── d2c/                         
│   ├── __init__.py
│   ├── ingredient2cell.py       
│   ├── herb2cell.py            
│   ├── formula2cell.py        
│   ├── ingreScan.py          
│   ├── herbScan.py          
│   ├── formulaScan.py      
│   ├── chemblScan.py            
│   └── ingreScanCell.py        
├── d4c/                       
│   ├── __init__.py
│   ├── TCMscore2cell.py        
│   └── TCMscoreScan.py        
├── pl/                        
│   ├── __init__.py
│   ├── TCMscoreScan_case_high.py  
│   ├── dotplot.py
│   └── umap.py
├── data/                        
├── utils/                      
│   ├── __init__.py
│   └──  data_check.py          
├── tests/                     
│   ├── bin                   
│   └── result               
├── setup.py                
├── pyproject.toml         
├── requirements.txt      
├── README.md            
└── docs/               
    ├── templates        
    └── output_format   


3. Templates
  The file templates used in the process of running sctcm are stored in the ./docs/templates directory.


4. Built-in data
  The built-in data required by sctcm is stored in the /sctcm/data directory.
