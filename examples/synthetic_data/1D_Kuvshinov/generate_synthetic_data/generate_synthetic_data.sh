./generate_MT_data setup_MT.config
./generate_Sq_G2LTF_data setup_Sq_G2LTF.config
./generate_Dst_Qn_data setup_Dst_Qn.config

python3 plot_synthetic_MT_data.py
python3 plot_synthetic_Sq_G2LTF_data.py
python3 plot_synthetic_Dst_Qn_data.py