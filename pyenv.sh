# pyenv environment 
python3.10 -m venv CSCN
source CSCN/bin/activate
pip install scanpy numpy pandas scikit-learn matplotlib pgmpy networkx jupyter
pip install fastapi uvicorn pytest
