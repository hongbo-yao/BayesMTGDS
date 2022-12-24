all: generate_Dst_Qn_data generate_Dst_Cn_data generate_Dst_G2LTF_data \
generate_Sq_Qn_data generate_Sq_Cn_data generate_Sq_G2LTF_data generate_MT_data

generate_Dst_Qn_data: generate_Dst_Qn_data.cpp ../src/utils.h ../src/forward_Qn.h
	g++ -o generate_Dst_Qn_data generate_Dst_Qn_data.cpp

generate_Dst_Cn_data: generate_Dst_Cn_data.cpp ../src/utils.h ../src/forward_Cn.h
	g++ -o generate_Dst_Cn_data generate_Dst_Cn_data.cpp

generate_Dst_G2LTF_data: generate_Dst_G2LTF_data.cpp ../src/utils.h ../src/forward_G2LTF.h
	g++ -o generate_Dst_G2LTF_data generate_Dst_G2LTF_data.cpp

generate_Sq_Qn_data: generate_Sq_Qn_data.cpp ../src/utils.h ../src/forward_Qn.h
	g++ -o generate_Sq_Qn_data generate_Sq_Qn_data.cpp

generate_Sq_Cn_data: generate_Sq_Cn_data.cpp ../src/utils.h ../src/forward_Cn.h
	g++ -o generate_Sq_Cn_data generate_Sq_Cn_data.cpp

generate_Sq_G2LTF_data: generate_Sq_G2LTF_data.cpp ../src/utils.h ../src/forward_G2LTF.h
	g++ -o generate_Sq_G2LTF_data generate_Sq_G2LTF_data.cpp

generate_MT_data: generate_MT_data.cpp ../src/utils.h ../src/forward_MT.h
	g++ -o generate_MT_data generate_MT_data.cpp

clean:
	rm generate_Dst_Qn_data generate_Dst_Cn_data generate_Dst_G2LTF_data \
	generate_Sq_Qn_data generate_Sq_Cn_data generate_Sq_G2LTF_data generate_MT_data