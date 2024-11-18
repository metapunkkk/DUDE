OPENQASM 2.0;
include "qelib1.inc";
qreg q[27];
creg c[27];
creg meas[27];
sx q[0];
rz(pi/4) q[0];
sx q[0];
sx q[1];
rz(0.61547971) q[1];
sx q[1];
sx q[2];
rz(pi/6) q[2];
sx q[2];
sx q[3];
rz(0.46364763) q[3];
sx q[3];
sx q[4];
rz(0.42053433) q[4];
sx q[4];
sx q[5];
rz(0.38759673) q[5];
sx q[5];
sx q[6];
rz(0.36136713) q[6];
sx q[6];
sx q[7];
rz(0.33983693) q[7];
sx q[7];
sx q[8];
rz(0.32175053) q[8];
sx q[8];
sx q[9];
rz(0.30627733) q[9];
sx q[9];
sx q[10];
rz(0.29284273) q[10];
sx q[10];
sx q[11];
rz(0.28103493) q[11];
sx q[11];
sx q[12];
rz(0.27054973) q[12];
sx q[12];
sx q[13];
rz(0.26115743) q[13];
sx q[13];
sx q[14];
rz(0.25268023) q[14];
sx q[14];
sx q[15];
rz(0.24497863) q[15];
sx q[15];
sx q[16];
rz(0.23794113) q[16];
sx q[16];
sx q[17];
rz(0.23147733) q[17];
sx q[17];
sx q[18];
rz(0.22551343) q[18];
sx q[18];
sx q[19];
rz(0.21998803) q[19];
sx q[19];
sx q[20];
rz(0.21484983) q[20];
sx q[20];
sx q[21];
rz(0.21005573) q[21];
sx q[21];
sx q[22];
rz(0.20556893) q[22];
sx q[22];
sx q[23];
rz(0.20135793) q[23];
sx q[23];
sx q[24];
rz(0.19739553) q[24];
sx q[24];
sx q[25];
rz(0.19365833) q[25];
sx q[25];
x q[26];
cx q[26],q[25];
sx q[25];
rz(0.19365833) q[25];
sx q[25];
cx q[25],q[24];
sx q[24];
rz(0.19739553) q[24];
sx q[24];
cx q[24],q[23];
sx q[23];
rz(0.20135793) q[23];
sx q[23];
cx q[23],q[22];
sx q[22];
rz(0.20556893) q[22];
sx q[22];
cx q[22],q[21];
sx q[21];
rz(0.21005573) q[21];
sx q[21];
cx q[21],q[20];
sx q[20];
rz(0.21484983) q[20];
sx q[20];
cx q[20],q[19];
sx q[19];
rz(0.21998803) q[19];
sx q[19];
cx q[19],q[18];
sx q[18];
rz(0.22551343) q[18];
sx q[18];
cx q[18],q[17];
sx q[17];
rz(0.23147733) q[17];
sx q[17];
cx q[17],q[16];
sx q[16];
rz(0.23794113) q[16];
sx q[16];
cx q[16],q[15];
sx q[15];
rz(0.24497863) q[15];
sx q[15];
cx q[15],q[14];
sx q[14];
rz(0.25268023) q[14];
sx q[14];
cx q[14],q[13];
sx q[13];
rz(0.26115743) q[13];
sx q[13];
cx q[13],q[12];
sx q[12];
rz(0.27054973) q[12];
sx q[12];
cx q[12],q[11];
sx q[11];
rz(0.28103493) q[11];
sx q[11];
cx q[11],q[10];
sx q[10];
rz(0.29284273) q[10];
sx q[10];
cx q[10],q[9];
cx q[25],q[26];
cx q[24],q[25];
cx q[23],q[24];
cx q[22],q[23];
cx q[21],q[22];
cx q[20],q[21];
cx q[19],q[20];
cx q[18],q[19];
cx q[17],q[18];
cx q[16],q[17];
cx q[15],q[16];
cx q[14],q[15];
cx q[13],q[14];
cx q[12],q[13];
cx q[11],q[12];
cx q[10],q[11];
sx q[9];
rz(0.30627733) q[9];
sx q[9];
cx q[9],q[8];
sx q[8];
rz(0.32175053) q[8];
sx q[8];
cx q[8],q[7];
sx q[7];
rz(0.33983693) q[7];
sx q[7];
cx q[7],q[6];
sx q[6];
rz(0.36136713) q[6];
sx q[6];
cx q[6],q[5];
sx q[5];
rz(0.38759673) q[5];
sx q[5];
cx q[5],q[4];
sx q[4];
rz(0.42053433) q[4];
sx q[4];
cx q[4],q[3];
sx q[3];
rz(0.46364763) q[3];
sx q[3];
cx q[3],q[2];
sx q[2];
rz(pi/6) q[2];
sx q[2];
cx q[2],q[1];
sx q[1];
rz(0.61547971) q[1];
sx q[1];
cx q[1],q[0];
sx q[0];
rz(pi/4) q[0];
sx q[0];
cx q[9],q[10];
cx q[8],q[9];
cx q[7],q[8];
cx q[6],q[7];
cx q[5],q[6];
cx q[4],q[5];
cx q[3],q[4];
cx q[2],q[3];
cx q[1],q[2];
cx q[0],q[1];