from flask import Flask, request, send_from_directory  
from flask_cors import CORS  
import os  
import time
from ase import Atoms
import ase.io
from ase.calculators.emt import EMT
import matgl
import torch
from matgl.ext.ase import M3GNetCalculator, Relaxer
from ase.calculators.singlepoint import SinglePointCalculator
from backend.utils.ase import Atoms2Structure

"""
Implemented Interfaces of M3/EGNet:
    - PES: e,f,s
    - Eform
    - BandGap
    - Relaxation
"""
  
app = Flask(__name__)  
CORS(app)  
  
UPLOAD_FOLDER = './uploads'
PROCESSED_FOLDER = './processed'
ALLOWED_EXTENSIONS = {'xyz', 'res', 'vasp', 'cif'}

os.system("rm -r {}".format(UPLOAD_FOLDER))
os.system("rm -r {}".format(PROCESSED_FOLDER))
os.makedirs(UPLOAD_FOLDER)
os.makedirs(PROCESSED_FOLDER)
  
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER  
app.config['PROCESSED_FOLDER'] = PROCESSED_FOLDER  
  
def allowed_file(filename):  
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS  
  
@app.route('/upload', methods=['POST'])  
def upload_file():  
    if request.method == 'POST':  
        file = request.files['file']  
        if file and allowed_file(file.filename):  
            filename = file.filename  
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))  
            result, res_file = process_file(filename)  
            os.remove(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            if result is not None:
                return {'status': 'success', 'filename': res_file, 'result': result}
            else:
                print("Calculation failed")
                return {'status': 'error', 'message': 'Calculation failed'}
        else:
            print("Invalid file type")
            return {'status': 'error', 'message': 'Invalid file type'}  
        
def process_file(filename):
    try:
        atoms = ase.io.read(os.path.join(UPLOAD_FOLDER, filename))
        struct = Atoms2Structure(atoms)
        pes_pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
        eform_model = matgl.load_model("M3GNet-MP-2018.6.1-Eform")
        bandgap_model = matgl.load_model("MEGNet-MP-2019.4.1-BandGap-mfi")
        calc = M3GNetCalculator(potential=pes_pot)
        atoms.set_calculator(calc)
        print("Computing energy, forces and stress...")
        energy = atoms.get_potential_energy().item()
        forces = atoms.get_forces().tolist()
        stress = atoms.get_stress().tolist()
        print("Computing eform...")
        eform = eform_model.predict_structure(struct).item()
        result = {
            "energy/eV": energy, ## eV
            # "forces": forces, ## eV/A
            # "stress": stress, ## TODO:eV/A^3 ???
            "eform/eV/atom": eform, ## eV/atom
        }
        print("Computing bandgap...")
        for i, method in ((0, "PBE"), (1, "GLLB-SC"), (2, "HSE"), (3, "SCAN")):
            graph_attrs = torch.tensor([i])
            bandgap = bandgap_model.predict_structure(
                structure=struct, state_feats=graph_attrs
            )
            result["BandGap_"+method+"/eV"] = bandgap.item() ## eV
        atoms.set_calculator(SinglePointCalculator(atoms, energy=energy, forces=forces, stress=stress))
        pure_name = filename.split('.')[0]
        ase.io.write(os.path.join(PROCESSED_FOLDER, pure_name+".xyz"), atoms)
    except:
        result = None
        pure_name = filename.split('.')[0]
    return result, pure_name+".xyz"
  
@app.route('/download/<path:filename>', methods=['GET'])  
def download_file(filename):
    return send_from_directory(directory=PROCESSED_FOLDER, path=filename, as_attachment=True)  

@app.route('/upload_relax', methods=['POST'])  
def upload_relax():
    print("Relaxing...")
    if request.method == 'POST':
        file = request.files['file']
        steps = int(request.form['steps'])
        fmax = float(request.form['fmax'])
        if file and allowed_file(file.filename):  
            filename = file.filename
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            try:
                pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
                relaxer = Relaxer(potential=pot)
                atoms = ase.io.read(os.path.join(UPLOAD_FOLDER, filename))
                struct = Atoms2Structure(atoms)
                pure_name = filename.split('.')[0]
                relax_results = relaxer.relax(
                    struct, 
                    steps=steps, 
                    fmax=fmax, 
                    traj_file=os.path.join(PROCESSED_FOLDER, "relax_{}.xyz".format(pure_name))
                )
                return {'status': 'success', 'filename': "relax_{}.xyz".format(pure_name)}
            except:
                print("Relaxation failed")
                return {'status': 'error', 'message': 'Relaxtion failed'}
        else:
            print("Invalid file type")
            return {'status': 'error', 'message': 'Invalid file type'}

if __name__ == '__main__':  
    app.run(debug=True)  
