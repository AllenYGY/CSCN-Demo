# conda environment file for CSCN
conda create -n CSCN python= 3.10
conda activate CSCN
conda install scanpy numpy pandas scikit-learn matplotlib pgmpy networkx -y
pip install fastapi uvicorn pytest
