// Project1.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <Windows.h>
#include <stdlib.h>
#include <gl/glut.h>
#include <math.h>
#include <string.h>
#include <vector>

#define judge 0.0001
#define e_judge 0.000000001
#define side_num 30

int node_Num_m[2];
int node_Num_n[2];
bool close_flag = false;
bool open_flag = false;
bool close_flag_n = false;
bool open_flag_n = false;
double c[3];
int window_size_x = 1000;
int window_size_y = window_size_x;
GLfloat orange[] = { 255.0 / 256.0, 153.0 / 256.0, 0.0 / 256.0, 0.9 };
GLfloat blue[] = { 0.0 / 256.0, 65.0 / 256.0, 255.0 / 256.0, 0.4 };
GLfloat white[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat blue2[] = { 102.0 / 256.0, 204.0 / 256.0, 255.0 / 256.0, 0.9 };
GLfloat blue_node[] = { 0.5, 0.5, 1.0, 1.0 };
bool up_flag = false;
bool down_flag = false;
int w_view;
int h_view;
int first_count = 1;
double wall_z = 0.5;
double wall_n = 5.0;
int num_count = 0;
int con_count = 0;
int tri_count = 0;
double damp_k = 1000.0;
double damp_k_normal = 20;

double dv = 3.0;
double node_Radius = 0.08;
double View_from[3] = { 0.0, 50.0, 50.0 };
double View_to[3] = { -50.0, 0.0, 0.0 };
double View_from2[3] = { 0.0, 13.0, 0.01 };
double View_to2[3] = { 0.0, -10.0, 0.0 };
double View_from3[3] = { 0.0, 13.0, 0.01 };
double View_to3[3] = { 0.0, -10.0, 0.0 };
bool MouseFlagRight = false;
bool MouseFlagLeft = false;
bool MouseFlagMiddle = false;
bool View_point_flag = false;

//追加変数
double last_mouse_x;
double last_mouse_y;
bool drag_mouse_r = false;
bool drag_mouse_l = false;
double ratio_x; //x,y軸を中心とする回転量
double ratio_y; 
double ratio_dis; //ズーム量
double camera_x = -30.0; //x,y軸を中心とする回転角度
double camera_y = -30.0;
double camera_dis = 15.0; //中心からのカメラの距離
				//

GLUnurbsObj *theNurb;
typedef struct {
	double x[3];
}position;
typedef struct {
	int number;
	int node_Num_w;
	int node_Num_h;
}face;
typedef struct {
	int t[3];
	int color;
	double normal[3];
	double A;
	double total;
}triangle_d;
typedef struct {
	bool torf;
	double len;
	//int color;
}edge_d;
typedef struct {
	int number;
	int edge_flag;			//ノード番号
	int none_flag;
	int dup_flag;
	int sin_flag;
	face face;
	position pos;		//位置
}node;
typedef struct {
	int number;
	int edge_flag;			//ノード番号
	int none_flag;
	int dup_flag;
	int sin_flag;
	face face;
	position pos;		//位置
	position del_pos;	//速度
	position acc;		//加速度
	double color_grad;
}node2;
typedef struct {
	int number;
	position pos;
	face face;
}point;
static node2 node_surface2[10000];
static node node_surface[1000][1000][2];
static edge_d edge[10000][10000];
static face face_info[3];
static point node_point[3][4]; //node_point[face Num][coordinate Num]
static triangle_d triangle_data[50000];
void mult_matrix3x1(double *c, double *a, double *b) {
	for (int i = 0; i < 3; i++) {
		c[i] = 0.0;
		for (int k = 0; k < 3; k++) {
			c[i] += a[3 * i + k] * b[k];
		}
	}
}
void gaiseki_9_3(double *a, double *b) {
	a[0] = 0.0;
	a[1] = -b[2];
	a[2] = b[1];
	a[3] = b[2];
	a[4] = 0.0;
	a[5] = -b[0];
	a[6] = -b[1];
	a[7] = b[0];
	a[8] = 0.0;
}
void sphere(double R, double precise, GLfloat sph_col[10]) {

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, sph_col);
	GLUquadricObj *sphere;

	sphere = gluNewQuadric();
	gluQuadricDrawStyle(sphere, GLU_FILL);

	gluSphere(sphere, R, precise, precise);
}
void View_control(double ratio_y) { //y軸（横
	// 現在の変換行列の右側に、今回の回転変換をかける
	glMatrixMode(GL_MODELVIEW);
	glRotatef(ratio_y, 0.0, 1.0, 0.0);
}
void View_control_up_down(double ratio_x) { //x軸（縦
	// 現在の変換行列を取得
	float  m[16];
	float  tx, ty, tz;
	glGetFloatv(GL_MODELVIEW_MATRIX, m);

	// 現在の変換行列の平行移動成分を記録
	tx = m[12];
	ty = m[13];
	tz = m[14];

	// 現在の変換行列の平行移動成分を０にする
	m[12] = 0.0f;
	m[13] = 0.0f;
	m[14] = 0.0f;

	// 変換行列を初期化
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// カメラの平行移動行列を設定
	glTranslatef(tx, ty, tz);

	// 右側に、今回の回転変換をかける
	glRotatef(ratio_x, 1.0, 0.0, 0.0);

	// さらに、右側に、もとの変換行列から平行移動成分をとり除いたものをかける
	glMultMatrixf(m);
}
void initiation() {

	int i = 0;
	int j = 0;
	int k = 0;
	int h = 0;
	int s = 0;
	int c = 0;
	int trirem1[3];
	int trirem2[3];
	int trirem3[3];
	int trirem4[3];
	int trirem5[3];
	int flag = 0;
	int tri_flag1;
	int tri_flag2;
	int tri_flag3;
	int tri_flag4;
	int tri_flag5;
	double tritemp_x;
	double tritemp_y;
	//ノードの間隔調節
	double natural_length = 1.0;
	//


	///三角にノードを置く
#if 1
	for (s = 0; s < 2; s++) {
		for (i = 0; i < side_num; i++) {
			for (j = 0; j <= i; j++) {
				if (i % 2 == 0) {
					node_surface[i][j][s].pos.x[0] = (double)((0.0) + (j - (i / 2)) * natural_length);
				}
				else {
					node_surface[i][j][s].pos.x[0] = (double)(((-i * 0.5) + j) * natural_length);
				}
				node_surface[i][j][s].pos.x[1] = (double)(i * sqrt(3) * 0.5);
				node_surface[i][j][s].pos.x[2] = (double)(0.0);
			}
		}
	}
#endif
	//

	//溶着面のノードの区別
#if 1
	for (s = 0; s < 2; s++) {
		for (i = 0; i < side_num; i++) {
			for (j = 0; j <= i; j++) {
				if (j == 0 || j == i || i == side_num - 1) {
					node_surface[i][j][s].edge_flag = 1;
					node_surface[i][j][s].none_flag = 0;
				}
				else {
					node_surface[i][j][s].edge_flag = 0;
					node_surface[i][j][s].none_flag = 1;
				}
			}
		}
	}
#endif
	//

	//ノードの番号付け
	for (s = 0; s < 2; s++) {
		for (i = 0; i < side_num; i++) {
			for (j = 0; j <= i; j++) {
				node_surface[i][j][s].number = num_count;
				num_count++;
			}
		}
	}
	printf("num_count = %d\n", num_count); //ノード数の表示										   
////////////////////////////////三角形のメッシュ切る
#if 1
	for (h = 0; h < num_count; h++) {
		for (s = 0; s < 2; s++) {
			tri_flag3 = 0;
			for (i = 0; i < side_num; i++) {
				for (j = 0; j <= i; j++) {
					if (node_surface[i][j][s].number == h) {
						trirem1[0] = i;
						trirem1[1] = j;
						trirem1[2] = s;
						tri_flag3 = 1;
					}
				}
			}
			tri_flag1 = 0;
			tri_flag2 = 0;
			tri_flag4 = 0;
			tri_flag5 = 0;
			for (i = 0; i < side_num; i++) {
				for (j = 0; j <= i; j++) {
					tritemp_x = node_surface[trirem1[0]][trirem1[1]][trirem1[2]].pos.x[0] - node_surface[i][j][s].pos.x[0];
					tritemp_y = node_surface[trirem1[0]][trirem1[1]][trirem1[2]].pos.x[1] - node_surface[i][j][s].pos.x[1];
					if (fabs(tritemp_x + natural_length * 0.5) < judge && fabs(tritemp_y + (sqrt(3) * natural_length * 0.5)) < judge) {
						trirem2[0] = i;
						trirem2[1] = j;
						trirem2[2] = s;
						tri_flag1 = 1;
					}
					if (fabs(tritemp_x - natural_length * 0.5) < judge && fabs(tritemp_y + (sqrt(3) * natural_length * 0.5)) < judge) {
						trirem3[0] = i;
						trirem3[1] = j;
						trirem3[2] = s;
						tri_flag2 = 1;
					}
					if (fabs(tritemp_x + natural_length) < judge && fabs(tritemp_y) < judge) {
						trirem4[0] = i;
						trirem4[1] = j;
						trirem4[2] = s;
						tri_flag4 = 1;
					}
					if (fabs(tritemp_x + natural_length * 0.5) < judge && fabs(tritemp_y + (sqrt(3) * natural_length * 0.5)) < judge) {
						trirem5[0] = i;
						trirem5[1] = j;
						trirem5[2] = s;
						tri_flag5 = 1;
					}
				}
			}
			if (tri_flag1 == 1 && tri_flag2 == 1 && tri_flag3 == 1) {
				if (tri_count == 0) {
					triangle_data[tri_count].t[0] = h;
					triangle_data[tri_count].t[1] = node_surface[trirem2[0]][trirem2[1]][trirem2[2]].number;
					triangle_data[tri_count].t[2] = node_surface[trirem3[0]][trirem3[1]][trirem3[2]].number;
					triangle_data[tri_count].color = 1;
					tri_count++;
				}
				else {
					flag = 0;
					for (i = 0; i < tri_count; i++) {
						if ((triangle_data[i].t[0] == h) && (triangle_data[i].t[1] == node_surface[trirem2[0]][trirem2[1]][trirem2[2]].number)
							&& (triangle_data[i].t[2] = node_surface[trirem3[0]][trirem3[1]][trirem3[2]].number)) {
							flag = 1;
						}
					}
					if (flag == 0) {
						triangle_data[tri_count].t[0] = h;
						triangle_data[tri_count].t[1] = node_surface[trirem2[0]][trirem2[1]][trirem2[2]].number;
						triangle_data[tri_count].t[2] = node_surface[trirem3[0]][trirem3[1]][trirem3[2]].number;
						triangle_data[tri_count].color = 1;
						tri_count++;
					}
					else {
						flag = 0;
					}
				}
			}
			if (tri_flag3 == 1 && tri_flag5 == 1 && tri_flag4 == 1) {
				if (tri_count == 0) {
					triangle_data[tri_count].t[0] = h;
					triangle_data[tri_count].t[1] = node_surface[trirem4[0]][trirem4[1]][trirem4[2]].number;
					triangle_data[tri_count].t[2] = node_surface[trirem5[0]][trirem5[1]][trirem5[2]].number;
					triangle_data[tri_count].color = 1;
					tri_count++;
				}
				else {
					flag = 0;
					for (i = 0; i < tri_count; i++) {
						if ((triangle_data[i].t[0] == h) && (triangle_data[i].t[1] == node_surface[trirem4[0]][trirem4[1]][trirem4[2]].number)
							&& (triangle_data[i].t[2] == node_surface[trirem5[0]][trirem5[1]][trirem5[2]].number)) {
							flag = 1;
						}
					}
					if (flag == 0) {
						triangle_data[tri_count].t[0] = h;
						triangle_data[tri_count].t[1] = node_surface[trirem4[0]][trirem4[1]][trirem4[2]].number;
						triangle_data[tri_count].t[2] = node_surface[trirem5[0]][trirem5[1]][trirem5[2]].number;
						triangle_data[tri_count].color = 1;
						tri_count++;
					}
					else {
						flag = 0;
					}
				}
			}
		}
	}
	//edge
#if 0
	for (s = 0; s < 2; s++) {
		for (i = 0; i < 13; i++) {
			for (j = 0; j <= i; j++) {
				if (s == 0 && node_surface[i][j][s].edge_flag == 1) {
					trirem1[0] = i;
					trirem1[1] = j;
					trirem1[2] = s;
				}
			}
		}
		for (k = 0; k < 3; k++) {
			node_surface[trirem1[0]][trirem1[1]][0].pos.x[k] = node_surface[trirem1[0]][trirem1[1]][1].pos.x[k];
		}
	}

#endif
#endif
	printf("%d\n", tri_count);

	//////////////////////////////////////node_surface2に同様に情報移転
#if 1
	for (s = 0; s < 2; s++) {
		for (i = 0; i < side_num; i++) {
			for (j = 0; j <= i; j++) {
				for (h = 0; h < 3; h++) {
					node_surface2[node_surface[i][j][s].number].pos.x[h] = node_surface[i][j][s].pos.x[h];
				}
			}
		}
	}
	for (s = 0; s < 2; s++) {
		for (i = 0; i < side_num; i++) {
			for (j = 0; j <= i; j++) {
				if (j == 0 || j == i || i == side_num - 1) {
					node_surface2[node_surface[i][j][s].number].edge_flag = 1;
					node_surface2[node_surface[i][j][s].number].none_flag = 0;
				}
				else {
					node_surface2[node_surface[i][j][s].number].edge_flag = 0;
					node_surface2[node_surface[i][j][s].number].none_flag = 1;
				}
			}
		}
	}
#endif	
	//////////////////////////////////////

	for (i = 0; i <= num_count - 1; i++) {
		for (j = 0; j <= num_count - 1; j++) {
			edge[i][j].torf = false;
		}
	}
	for (i = 0; i < tri_count; i++) {
		edge[triangle_data[i].t[0]][triangle_data[i].t[1]].torf = true;
		edge[triangle_data[i].t[1]][triangle_data[i].t[0]].torf = true;
		edge[triangle_data[i].t[1]][triangle_data[i].t[2]].torf = true;
		edge[triangle_data[i].t[2]][triangle_data[i].t[1]].torf = true;
		edge[triangle_data[i].t[2]][triangle_data[i].t[0]].torf = true;
		edge[triangle_data[i].t[0]][triangle_data[i].t[2]].torf = true;
		edge[triangle_data[i].t[0]][triangle_data[i].t[1]].len = sqrt(
			pow((node_surface2[triangle_data[i].t[0]].pos.x[0] - node_surface2[triangle_data[i].t[1]].pos.x[0]), 2.0) +
			pow((node_surface2[triangle_data[i].t[0]].pos.x[1] - node_surface2[triangle_data[i].t[1]].pos.x[1]), 2.0) +
			pow((node_surface2[triangle_data[i].t[0]].pos.x[2] - node_surface2[triangle_data[i].t[1]].pos.x[2]), 2.0));
		edge[triangle_data[i].t[1]][triangle_data[i].t[0]].len = sqrt(
			pow((node_surface2[triangle_data[i].t[1]].pos.x[0] - node_surface2[triangle_data[i].t[0]].pos.x[0]), 2.0) +
			pow((node_surface2[triangle_data[i].t[1]].pos.x[1] - node_surface2[triangle_data[i].t[0]].pos.x[1]), 2.0) +
			pow((node_surface2[triangle_data[i].t[1]].pos.x[2] - node_surface2[triangle_data[i].t[0]].pos.x[2]), 2.0));
		edge[triangle_data[i].t[1]][triangle_data[i].t[2]].len = sqrt(
			pow((node_surface2[triangle_data[i].t[1]].pos.x[0] - node_surface2[triangle_data[i].t[2]].pos.x[0]), 2.0) +
			pow((node_surface2[triangle_data[i].t[1]].pos.x[1] - node_surface2[triangle_data[i].t[2]].pos.x[1]), 2.0) +
			pow((node_surface2[triangle_data[i].t[1]].pos.x[2] - node_surface2[triangle_data[i].t[2]].pos.x[2]), 2.0));
		edge[triangle_data[i].t[2]][triangle_data[i].t[1]].len = sqrt(
			pow((node_surface2[triangle_data[i].t[2]].pos.x[0] - node_surface2[triangle_data[i].t[1]].pos.x[0]), 2.0) +
			pow((node_surface2[triangle_data[i].t[2]].pos.x[1] - node_surface2[triangle_data[i].t[1]].pos.x[1]), 2.0) +
			pow((node_surface2[triangle_data[i].t[2]].pos.x[2] - node_surface2[triangle_data[i].t[1]].pos.x[2]), 2.0));
		edge[triangle_data[i].t[2]][triangle_data[i].t[0]].len = sqrt(
			pow((node_surface2[triangle_data[i].t[2]].pos.x[0] - node_surface2[triangle_data[i].t[0]].pos.x[0]), 2.0) +
			pow((node_surface2[triangle_data[i].t[2]].pos.x[1] - node_surface2[triangle_data[i].t[0]].pos.x[1]), 2.0) +
			pow((node_surface2[triangle_data[i].t[2]].pos.x[2] - node_surface2[triangle_data[i].t[0]].pos.x[2]), 2.0));
		edge[triangle_data[i].t[0]][triangle_data[i].t[2]].len = sqrt(
			pow((node_surface2[triangle_data[i].t[0]].pos.x[0] - node_surface2[triangle_data[i].t[2]].pos.x[0]), 2.0) +
			pow((node_surface2[triangle_data[i].t[0]].pos.x[1] - node_surface2[triangle_data[i].t[2]].pos.x[1]), 2.0) +
			pow((node_surface2[triangle_data[i].t[0]].pos.x[2] - node_surface2[triangle_data[i].t[2]].pos.x[2]), 2.0));
	}
}
void node_simulation(int view_con) {
	if (first_count == 1) {
		initiation();
		first_count--;
	}
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	int h = 0;
	int s = 0;
	int E = i + num_count * 0.5;
	double min = 100000.0;
	double max = -10.0;
	double mass = 1000;
	double kizami = 0.01; //時間t(0.01秒ごと+)
	double total_force[3] = { 0.0, 0.0, 0.0 };
	double temp_len;
	double normal_force[3];
	double normal_temp[9];
	double normal_temp3[3];
	double pressure;

	//全ての点の各方向の加速度を0に設定
	for (i = 0; i <= num_count - 1; i++) {
		for (j = 0; j < 3; j++) {
			node_surface2[i].acc.x[j] = 0.0;
		}
	}


	for (i = 0; i <= num_count - 1; i++) {
		node_surface2[i].color_grad = 0.0;
		for (j = 0; j < num_count; j++) {
			if (edge[i][j].torf == 1) {
				temp_len = sqrt(pow((node_surface2[i].pos.x[0] - node_surface2[j].pos.x[0]), 2.0) + pow((node_surface2[i].pos.x[1] - node_surface2[j].pos.x[1]), 2.0) + pow((node_surface2[i].pos.x[2] - node_surface2[j].pos.x[2]), 2.0));
				if (temp_len > edge[i][j].len) {
					for (k = 0; k < 3; k++) {
						node_surface2[i].color_grad += fabs(-1000.0 * damp_k * (node_surface2[i].pos.x[k] - node_surface2[j].pos.x[k]) * (temp_len - edge[i][j].len) / temp_len / mass);
						node_surface2[i].acc.x[k] += -1000.0 * damp_k * (node_surface2[i].pos.x[k] - node_surface2[j].pos.x[k]) * (temp_len - edge[i][j].len) / temp_len / mass;
					}
				}
			}
		}
	}

	if (open_flag == true) {
		wall_z += 0.01;
	}
	if (close_flag == true) {
		wall_z -= 0.01;
	}
	for (i = 0; i < num_count; i++) {
		for (k = 0; k < 3; k++) {
			node_surface2[i].acc.x[k] += -dv * node_surface2[i].del_pos.x[k];
		}
	}
#if 1

	for (i = 0; i < tri_count; i++) {
		if (i < tri_count * 0.5) {
			if (triangle_data[i].color == 1) {
				for (j = 0; j < 3; j++) {
					normal_temp3[j] = node_surface2[triangle_data[i].t[1]].pos.x[j] - node_surface2[triangle_data[i].t[0]].pos.x[j];
				}
				gaiseki_9_3(normal_temp, normal_temp3);
				for (j = 0; j < 3; j++) {
					normal_temp3[j] = node_surface2[triangle_data[i].t[2]].pos.x[j] - node_surface2[triangle_data[i].t[0]].pos.x[j];
				}
				mult_matrix3x1(normal_force, normal_temp, normal_temp3);
			}
			else if (triangle_data[i].color == 2) {
				for (j = 0; j < 3; j++) {
					normal_temp3[j] = node_surface2[triangle_data[i].t[2]].pos.x[j] - node_surface2[triangle_data[i].t[0]].pos.x[j];
				}
				gaiseki_9_3(normal_temp, normal_temp3);
				for (j = 0; j < 3; j++) {
					normal_temp3[j] = node_surface2[triangle_data[i].t[1]].pos.x[j] - node_surface2[triangle_data[i].t[0]].pos.x[j];
				}
				mult_matrix3x1(normal_force, normal_temp, normal_temp3);
			}
		}
		if (i >= tri_count * 0.5) {
			if (triangle_data[i].color == 1) {
				for (j = 0; j < 3; j++) {
					normal_temp3[j] = node_surface2[triangle_data[i].t[2]].pos.x[j] - node_surface2[triangle_data[i].t[0]].pos.x[j];
				}
				gaiseki_9_3(normal_temp, normal_temp3);
				for (j = 0; j < 3; j++) {
					normal_temp3[j] = node_surface2[triangle_data[i].t[1]].pos.x[j] - node_surface2[triangle_data[i].t[0]].pos.x[j];
				}
				mult_matrix3x1(normal_force, normal_temp, normal_temp3);
			}
			else if (triangle_data[i].color == 2) {
				for (j = 0; j < 3; j++) {
					normal_temp3[j] = node_surface2[triangle_data[i].t[1]].pos.x[j] - node_surface2[triangle_data[i].t[0]].pos.x[j];
				}
				gaiseki_9_3(normal_temp, normal_temp3);
				for (j = 0; j < 3; j++) {
					normal_temp3[j] = node_surface2[triangle_data[i].t[2]].pos.x[j] - node_surface2[triangle_data[i].t[0]].pos.x[j];
				}
				mult_matrix3x1(normal_force, normal_temp, normal_temp3);
			}
		}
		for (j = 0; j < 3; j++) {
			triangle_data[i].normal[j] = normal_force[j] / sqrt(pow(normal_force[0], 2.0) + pow(normal_force[1], 2.0) + pow(normal_force[2], 2.0));// / normal_force[j];
			for (k = 0; k < 3; k++) {
				node_surface2[triangle_data[i].t[j]].acc.x[k] += 0.1 * damp_k_normal * normal_force[k];
			}
		}
	}
#endif
#if 1
	//質点の位置、速度の更新(node_surface2_none_flag)
	for (i = 0; i < num_count; i++) {
		if (node_surface2[i].none_flag == 1) {
			for (k = 0; k < 3; k++) {
				node_surface2[i].del_pos.x[k] = node_surface2[i].del_pos.x[k] + node_surface2[i].acc.x[k] * kizami; //速度の更新
				node_surface2[i].pos.x[k] = node_surface2[i].pos.x[k] + node_surface2[i].del_pos.x[k] * kizami + 1.0 / 2.0 * node_surface2[i].acc.x[k] * kizami * kizami; //質点の位置の更新
			}
		}
	}

	//質点の位置、速度の更新(node_surface2_edge_flag)
	for (i = 0; i < num_count * 0.5; i++) {
		if (node_surface2[i].edge_flag == 1) {
			for (k = 0; k < 3; k++) {
				int j = i + num_count * 0.5;
				node_surface2[i].del_pos.x[k] = node_surface2[i].del_pos.x[k] + (node_surface2[i].acc.x[k] * 0.5 + node_surface2[j].acc.x[k] * 0.5) * kizami;
				node_surface2[j].del_pos.x[k] = node_surface2[j].del_pos.x[k] + (node_surface2[i].acc.x[k] * 0.5 + node_surface2[j].acc.x[k] * 0.5) * kizami;
				node_surface2[i].pos.x[k] = node_surface2[i].pos.x[k] + node_surface2[i].del_pos.x[k] * kizami + 1.0 / 4.0 * node_surface2[i].acc.x[k] * kizami * kizami + node_surface2[j].del_pos.x[k] * kizami + 1.0 / 4.0 * node_surface2[j].acc.x[k] * kizami * kizami;
				node_surface2[j].pos.x[k] = node_surface2[j].pos.x[k] + node_surface2[i].del_pos.x[k] * kizami + 1.0 / 4.0 * node_surface2[i].acc.x[k] * kizami * kizami + node_surface2[j].del_pos.x[k] * kizami + 1.0 / 4.0 * node_surface2[j].acc.x[k] * kizami * kizami;
			}
		}
	}

#endif
	GLfloat changing[] = { 0.5, 0.5, 1.0, 1.0 }; // Blue

	for (i = 0; i < num_count; i++) {
		if (node_surface2[i].color_grad > max) {
			max = node_surface2[i].color_grad;
		}
		if (node_surface2[i].color_grad < min) {
			min = node_surface2[i].color_grad;
		}
	}
	for (i = 0; i < num_count; i++) {
		glPushMatrix();
		glCullFace(GL_BACK);
		if (node_surface2[i].color_grad < (max - min) / 2.0) {
			changing[0] = (node_surface2[i].color_grad - min) / ((max + min) / 2.0 - min);
			changing[1] = 0.1;
			changing[2] = ((max + min) / 2.0 - node_surface2[i].color_grad) / ((max + min) / 2.0 - min);
		}
		else {
			changing[0] = (max - node_surface2[i].color_grad) / ((max + min) / 2.0 - min);
			changing[1] = (node_surface2[i].color_grad - (max + min) / 2.0) / ((max + min) / 2.0 - min);
			changing[2] = 0.1;
		}
		glTranslated((GLdouble)node_surface2[i].pos.x[0] - 1, (GLdouble)node_surface2[i].pos.x[2], (GLdouble)node_surface2[i].pos.x[1] - side_num * 0.5);
		if (view_con == 1) glutSolidSphere(node_Radius, 10, 10);
		else if (view_con == 2) {
			glMaterialfv(GL_FRONT, GL_DIFFUSE, changing);
			glutSolidCube(node_Radius * 4.0);
		}
		glPopMatrix();
	}
	glPushMatrix();
	//glTranslated( - rec_x / 2.0,0.0,  - rec_y / 2.0);
	for (i = 0; i < tri_count; i++) {
		if (i < tri_count * 0.5) {
			if (triangle_data[i].color == 1) {
				glCullFace(GL_FRONT);
			}
			else {
				glCullFace(GL_BACK);
			}
		}
		else if (i >= tri_count * 0.5) {
			if (triangle_data[i].color == 1) {
				glCullFace(GL_BACK);
			}
			else {
				glCullFace(GL_FRONT);
			}
		}
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		if (view_con == 2) {
			glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
		}
		else {
			glMaterialfv(GL_FRONT, GL_DIFFUSE, blue2);
		}
		//物体を描画(三角形)
		glBegin(GL_TRIANGLES);
		//printf("%f %f %f %f\n", blue[0], blue[1], blue[2], blue[3]);
		if (1) {
			glNormal3d(triangle_data[i].normal[0], triangle_data[i].normal[2], triangle_data[i].normal[1]);
			glVertex3d(node_surface2[triangle_data[i].t[0]].pos.x[0] - 1, node_surface2[triangle_data[i].t[0]].pos.x[2], node_surface2[triangle_data[i].t[0]].pos.x[1] - side_num * 0.5);
			glVertex3d(node_surface2[triangle_data[i].t[1]].pos.x[0] - 1, node_surface2[triangle_data[i].t[1]].pos.x[2], node_surface2[triangle_data[i].t[1]].pos.x[1] - side_num * 0.5);
			glVertex3d(node_surface2[triangle_data[i].t[2]].pos.x[0] - 1, node_surface2[triangle_data[i].t[2]].pos.x[2], node_surface2[triangle_data[i].t[2]].pos.x[1] - side_num * 0.5);
		}
		glEnd();
	}
	glPopMatrix();
}
void idle(void)
{
	glutPostRedisplay();
}
void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	//
	glTranslatef(0.0, 0.0, -camera_dis);
	glRotatef(-camera_x, 1.0, 0.0, 0.0);
	glRotatef(-camera_y, 0.0, 1.0, 0.0);
	//
	gluLookAt(View_from[0], View_from[1], View_from[2],
		View_to[0], View_to[1], View_to[2],
		0.0, 1.0, 0.0);
	glViewport(0, 0, w_view * 2.0 / 3.0, h_view);
	glPushMatrix();
	node_simulation(1);
	glPopMatrix();

	glutSwapBuffers();

}
//新マウス操作(ドラッグの検出)
void mouse(int button, int state, int x, int y) {
	//右ボタンのドラッグ検出
	if ((button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN)) {
		drag_mouse_l = 1;
	}
	else if ((button == GLUT_LEFT_BUTTON) && (state == GLUT_UP)) {
		drag_mouse_l = 0;
	}

	//左ボタンのドラッグ検出
	if ((button == GLUT_RIGHT_BUTTON) && (state == GLUT_DOWN)) {
		drag_mouse_r = 1;
	}
	else if ((button == GLUT_RIGHT_BUTTON) && (state == GLUT_UP)) {
		drag_mouse_r = 0;
	}

	//現在のマウス座標の記録
	last_mouse_x = x;
	last_mouse_y = y;
	/*
	switch (button) { //ドラッグの検知
	case GLUT_LEFT_BUTTON:
		if (state == GLUT_DOWN) {
			drag_mouse_l = 1;
			last_mouse_x = x; //マウスの最終位置(click時)の記録
			last_mouse_y = y;
		}
		else drag_mouse_l = 0;
		break;
	case GLUT_RIGHT_BUTTON:
		if (state == GLUT_DOWN) {
			drag_mouse_r = 1;
			last_mouse_z = x;
		}
		else drag_mouse_r = 0;
		break;
	default:
		break;
	}
	*/
}
//ドラッグ時の挙動
void motion(int x, int y) {
	if (drag_mouse_l == 1) {
		//y軸中心回転（横
		ratio_y -= (x - last_mouse_x) * 1.0;
		if (camera_x < -90.0) camera_x = -90.0;
		//x軸中心回転（縦
		ratio_x -= (y - last_mouse_y) * 1.0;
		if (camera_y < 0.0) camera_y += 360.0;
		else if (camera_y > 360.0) camera_y -= 360.0;
	}
	//z軸
	else if (drag_mouse_r == 1) {
		//マウスの横移動に応じてカメラの距離を移動
		ratio_dis = (x - last_mouse_x) * 1.0;
		//マウスの横移動に応じて距離を移動
		camera_dis += ratio_dis * 0.01;
		if (camera_dis < 5.0) camera_dis = 5.0;
		//現在の変換行列（カメラの向き）を取得
		float m[16];
		glGetFloatv(GL_MODELVIEW_MATRIX, m);

		// 変換行列を初期化して、カメラ移動分の平行移動行列を設定
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslatef(0.0, 0.0, -ratio_dis);

		// 右からこれまでの変換行列をかける
		glMultMatrixf(m);

	}
	//カメラ操作
	if (drag_mouse_l) {
		View_control(ratio_y);
		View_control_up_down(ratio_x);
	}
	//現在のマウス位置の記録
	last_mouse_x = x;
	last_mouse_y = y;
}

void resize(int w, int h)
{
	w_view = w;
	h_view = h;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(50.0, (double)w / (double)h * 2.0 / 3.0, 1.0, 1000.0);

	glMatrixMode(GL_MODELVIEW);

}
void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 'r':
		close_flag = true;
		break;
	case 'e':
		close_flag = false;
		break;
	case 'w':
		open_flag = true;
		break;
	case 'q':
		open_flag = false;
		break;
	case 'y':
		close_flag_n = true;
		break;
	case 'u':
		close_flag_n = false;
		break;
	case 'i':
		open_flag_n = true;
		break;
	case 'o':
		open_flag_n = false;
		break;
	}
}
void init() {

	glClearColor(1.0, 1.0, 1.0, 1.0);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	GLfloat direction[] = { 0.0, 1.0, 0.0 };
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);//アルファの設定
	glEnable(GL_BLEND);//アルファのブレンド有効
	glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, direction);
	glLightfv(GL_LIGHTING, GL_SPOT_DIRECTION, direction);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	glDisable(GL_LIGHT4);

}

int main(int argc, char *argv[])
{

	glutInit(&argc, argv);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(window_size_x, window_size_y);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow(argv[0]);
	glutInitWindowPosition(0, 0);
	glutDisplayFunc(display);
	glutReshapeFunc(resize);
	glutIdleFunc(idle);
	glutMouseFunc(mouse);
	glutMotionFunc(motion); //追加(ドラッグ時の挙動)
	glutKeyboardFunc(keyboard);
	init();
	glutMainLoop();

	return 0;

}