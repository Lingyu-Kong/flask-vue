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
from backend.utils.ase import Atoms2Structure, Structure2Atoms
import shutil  
import schedule  
import time  
from threading import Thread  

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

def clean_up():  
    folders = ['UNPROCESSED', 'PROCESSED']  
      
    current_time = time.time()  
  
    for folder in folders:  
        for filename in os.listdir(folder):  
            file_path = os.path.join(folder, filename)  
            try:  
                # 计算文件的存放时间（以秒为单位）  
                file_age = current_time - os.path.getmtime(file_path)  
  
                # 如果文件已存在超过24小时（24*60*60秒），则删除文件  
                if file_age > 24 * 60 * 60:  
                    if os.path.isfile(file_path) or os.path.islink(file_path):  
                        os.unlink(file_path)  
                    elif os.path.isdir(file_path):  
                        shutil.rmtree(file_path)  
            except Exception as e:  
                print(f'Failed to delete {file_path}. Reason: {e}')  
  
def run_schedule():  
    while True:  
        schedule.run_pending()  
        time.sleep(1)  
  
# 在此处设置定期任务的时间间隔，例如：每小时执行一次  
schedule.every(1).hours.do(clean_up)  
  
# 使用一个单独的线程运行定时任务  
t = Thread(target=run_schedule)  
t.start()  
  
def allowed_file(filename):  
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS  
  
@app.route('/pes_calc', methods=['POST'])
def single_point_calculation(): 
    if request.method == 'POST':  
        file = request.files['file']  
        if file and allowed_file(file.filename):  
            filename = file.filename  
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))  
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
                    # "energy/eV": energy, ## eV
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
            res_file = pure_name+".xyz"
            os.remove(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            if result is not None:
                return {'status': 'success', 'res_file': res_file, 'result': result}
            else:
                print("Calculation failed")
                return {'status': 'error', 'message': 'Calculation failed'}
        else:
            print("Invalid file type")
            return {'status': 'error', 'message': 'Invalid file type'}   

@app.route('/relax', methods=['POST'])  
def relax():
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
                    traj_file=os.path.join(PROCESSED_FOLDER, "{}_traj.xyz".format(pure_name))
                )
                print("Relaxation finished")
                relaxed_struct = relax_results["final_structure"]
                relaxed_atoms = Structure2Atoms(relaxed_struct)
                pes_pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
                eform_model = matgl.load_model("M3GNet-MP-2018.6.1-Eform")
                bandgap_model = matgl.load_model("MEGNet-MP-2019.4.1-BandGap-mfi")
                calc = M3GNetCalculator(potential=pes_pot)
                relaxed_atoms.set_calculator(calc)
                print("Computing energy, forces and stress...")
                energy = relaxed_atoms.get_potential_energy().item()
                forces = relaxed_atoms.get_forces().tolist()
                stress = relaxed_atoms.get_stress().tolist()
                print("Computing eform...")
                eform = eform_model.predict_structure(relaxed_struct).item()
                result = {
                    # "energy/eV": energy, ## eV
                    # "forces": forces, ## eV/A
                    # "stress": stress, ## TODO:eV/A^3 ???
                    "eform/eV/atom": eform, ## eV/atom
                }
                print("Computing bandgap...")
                for i, method in ((0, "PBE"), (1, "GLLB-SC"), (2, "HSE"), (3, "SCAN")):
                    graph_attrs = torch.tensor([i])
                    bandgap = bandgap_model.predict_structure(
                        structure=relaxed_struct, state_feats=graph_attrs
                    )
                    result["BandGap_"+method+"/eV"] = bandgap.item() ## eV
                relaxed_atoms.set_calculator(SinglePointCalculator(relaxed_atoms, energy=energy, forces=forces, stress=stress))
                pure_name = filename.split('.')[0]
                ase.io.write(os.path.join(PROCESSED_FOLDER, "{}_relaxed.xyz".format(pure_name)), atoms)
                os.remove(os.path.join(app.config['UPLOAD_FOLDER'], filename))
                return {
                    'status': 'success', 
                    'relaxed_file': "{}_relaxed.xyz".format(pure_name),
                    'traj_file': "{}_traj.xyz".format(pure_name),
                    'result': result
                }
            except:
                print("Relaxation failed")
                return {'status': 'error', 'message': 'Relaxtion failed'}
        else:
            print("Invalid file type")
            return {'status': 'error', 'message': 'Invalid file type'}
        
@app.route('/download/<path:filename>', methods=['GET'])  
def download_file(filename):
    return send_from_directory(directory=PROCESSED_FOLDER, path=filename, as_attachment=True) 

if __name__ == '__main__':  
    app.run(debug=True)  