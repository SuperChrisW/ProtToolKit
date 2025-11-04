#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import csv
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.DSSP import DSSP

class ProteinSelect(Select):
    """用于选择要保留的原子的类"""
    def accept_residue(self, residue):
        # 排除水分子和非标准氨基酸
        hetero, resseq, icode = residue.get_id()
        if hetero == "W":  # 水分子
            return 0
        if residue.get_resname() in ["HOH", "WAT"]:  # 水分子
            return 0
        if hetero != " ":  # 其他非标准分子（如配体等）
            return 0
        return 1

def standardize_pdb(pdb_file, output_dir, missing_data):
    """标准化单个PDB文件并检查缺失原子"""
    pdb_id = os.path.basename(pdb_file).split('.')[0]
    parser = PDBParser(QUIET=True)
    
    try:
        structure = parser.get_structure(pdb_id, pdb_file)
        
        # 检查缺失原子
        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                for residue in chain:
                    if residue.get_id()[0] == " ":  # 只检查标准氨基酸
                        res_name = residue.get_resname()
                        res_id = residue.get_id()[1]
                        
                        # 检查主链原子
                        for atom_name in ["N", "CA", "C", "O"]:
                            if atom_name not in residue:
                                missing_data.append([
                                    pdb_id,
                                    chain_id,
                                    res_name,
                                    res_id,
                                    atom_name
                                ])
        
        # 标准化链标识符
        # 如果链ID为空或非标准，则重新编号
        chain_mapping = {}
        new_chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 
                      'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
        chain_index = 0
        
        for model in structure:
            for chain in model:
                old_id = chain.get_id()
                if old_id == " " or not old_id.isalpha():
                    new_id = new_chains[chain_index]
                    chain_index += 1
                    chain.id = new_id
                    chain_mapping[old_id] = new_id
        
        # 保存标准化后的文件
        output_file = os.path.join(output_dir, f"{pdb_id}_standardized.pdb")
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_file, ProteinSelect())
        
        return True
    except Exception as e:
        print(f"处理 {pdb_id} 时出错: {str(e)}")
        return False

def process_all_pdbs(input_dir, output_dir):
    """处理目录中的所有PDB文件"""
    # 确保输出目录存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 获取所有PDB文件
    pdb_files = glob.glob(os.path.join(input_dir, "*.pdb")) + glob.glob(os.path.join(input_dir, "*.ent"))
    
    if not pdb_files:
        print(f"警告: 在 {input_dir} 中没有找到PDB文件")
        return
    
    # 用于存储缺失原子信息的列表
    missing_data = []
    
    # 处理每个PDB文件
    success_count = 0
    for pdb_file in pdb_files:
        print(f"正在处理: {os.path.basename(pdb_file)}")
        if standardize_pdb(pdb_file, output_dir, missing_data):
            success_count += 1
    
    # 将缺失原子信息写入CSV文件
    csv_file = os.path.join(output_dir, "missing_atoms.csv")
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["PDB_ID", "链", "残基名", "残基编号", "缺失原子"])
        writer.writerows(missing_data)
    
    print(f"\n处理完成!")
    print(f"成功处理: {success_count}/{len(pdb_files)} 个PDB文件")
    print(f"发现 {len(missing_data)} 个缺失原子")
    print(f"缺失原子报告已保存到: {csv_file}")
    print(f"标准化后的PDB文件已保存到: {output_dir}")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="标准化PDB文件并检查缺失原子")
    parser.add_argument("input_dir", help="包含PDB文件的输入目录")
    parser.add_argument("--output_dir", "-o", default="standardized_pdbs", 
                        help="保存标准化PDB文件的输出目录")
    
    args = parser.parse_args()
    
    process_all_pdbs(args.input_dir, args.output_dir)
