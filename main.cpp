#include <cstdio>
#include <array>
#include <vector>
#include <random>
#include "AVLTreeSum3.h"
#include "KMCLogger.h"
#include <fstream>
#include <iostream>
#include <chrono>

inline static const double kbT = 0.5;  //0.034から0.394?//
inline static const double D = 0.347;
inline static const double C = (8.9 - 8 * D) / sqrt(8);
inline static const double material_density = 1.0;
inline static const double material_density_2 = 1.0;

struct MoveTarget {
	int size_id_org;      //移動前のサイトを示すID//
	int size_id_dest;     //移動後のサイトを示すID//
	double ratio;         //移動頻度(単位時間当たりの発生確率)//
};

//あるサイトから見て隣(移動先)となるサイト(パス)の数の上限//
//ユーザーの設定する格子に依存して適切な値を設定せよ//

//constexpr size_t NUM_NEIGHBOR_SITES = 4; 
struct EventAtom {//同一の粒子でも複数の移動先を持つので、それらを一つの情報に纏める為のクラス//
	int atom_id;            //粒子固有のID
	int currest_site_id;    //粒子が現在位置しているサイト//
	std::vector<MoveTarget> paths;
	//std::array<MoveTarget, NUM_NEIGHBOR_SITES> paths;
	//int num_effectivepaths = 0;    //有効な移動先の数//

	double Ratio() {
		double ratio = 0.0;
		for (const auto& a : paths) {
			ratio += a.ratio;
		}
		return ratio;
	}
};

using KMCEventList = AVLTreeSum3<int, double, EventAtom*>;

/*********************************************
* 基本ルール
* - 粒子は決められたサイト(格子点)の上を移動するものとする
* - サイトは場所ごとに一意に決められた固有のIDを持つ
*   - 逆にIDが分かればサイトの座標などが分かるものとする
* 
***********************************************/


struct SiteInfo {
	int exist_atom_id;  //サイトに存在する粒子のID//存在しないときはATOM_NONEとなる
	static const int ATOM_NONE = -1;
};

/**
* 初期の粒子配置を決める
* 計算したい材料の構造に合わせて適宜書き換える
* この例ではsitesで与えられているサイトに、ある確率に従って粒子を配置する
**/
void InitializeSite(std::vector<SiteInfo>& sites, std::vector<EventAtom>& atoms, int lattice_x, int lattice_y, int lattice_z) {
	unsigned int seed = 123456789;
	std::mt19937 mt(seed);
	std::uniform_real_distribution<> dist(0.0, 1.0);

	int num_atoms = 0;
	int site_id = 0;
	int N = sites.size();
	int M = lattice_x * lattice_y * 20  ; // 下から1000個分の粒子

	// x方向を5等分
	int region_width_x = lattice_x / 10;
	// z方向を5等分 (z > 10の領域のみ)
	int region_height_z = (lattice_z - 10) / 10;

	// 各領域の粒子配置フラグ (5x5の2次元配列)
	std::vector<std::vector<bool>> region_flags(10, std::vector<bool>(10, true));

	// ここで各領域のフラグを設定します
	//1段目
	//region_flags[0][0] = false;
	//region_flags[1][0] = false;
	//region_flags[2][0] = false;
	//region_flags[3][0] = false;
	//region_flags[4][0] = false;
	//region_flags[5][0] = false;

	//2段目
	//region_flags[0][1] = false;
	//region_flags[1][1] = false;
	//region_flags[2][1] = false;
	//region_flags[3][1] = false;
	//region_flags[4][1] = false;

	//3段目
	//region_flags[0][2] = false;
	//region_flags[1][2] = false;
    //region_flags[2][2] = false;
	//region_flags[3][2] = false;
	//region_flags[4][2] = false;
	//region_flags[5][2] = false;
	//region_flags[6][2] = false;
	//region_flags[7][2] = false;
	//region_flags[8][2] = false;
	//region_flags[9][2] = false;
	//4段目
	//region_flags[0][3] = false;
	//region_flags[1][3] = false;
	//region_flags[2][3] = false;
	//region_flags[3][3] = false;
	//region_flags[4][3] = false;
	//region_flags[5][3] = false;
	//region_flags[6][3] = false;
	//region_flags[7][3] = false;
	//region_flags[8][3] = false;
	//region_flags[9][3] = false;

	//5段目
	//region_flags[0][4] = false;
	//region_flags[1][4] = false;
	region_flags[2][4] = false;
	region_flags[3][4] = false;
	//region_flags[4][4] = false;
	region_flags[5][4] = false;
	region_flags[6][4] = false;
	//region_flags[7][4] = false;
	//region_flags[9][4] = false;

	//6段目
	//region_flags[0][5] = false;
	//region_flags[1][5] = false;
	region_flags[2][5] = false;
	region_flags[3][5] = false;
	//region_flags[4][5] = false;
	region_flags[5][5] = false;
	region_flags[6][5] = false;
	//region_flags[7][5] = false;
	//region_flags[9][5] = false;
	
	//7段目
	//region_flags[0][6] = false;
	//region_flags[1][6] = false;
	//region_flags[2][6] = false;
	//region_flags[3][6] = false;
	//region_flags[4][6] = false;
	//region_flags[5][6] = false;

	//region_flags[9][6] = false;

	//8段目
	/*
	region_flags[0][7] = false;
	region_flags[1][7] = false;
	region_flags[2][7] = false;
	region_flags[3][7] = false;
	region_flags[4][7] = false;
	region_flags[5][7] = false;
	region_flags[6][7] = false;
	region_flags[7][7] = false;
	region_flags[8][7] = false;
	region_flags[9][7] = false;
	

	//9段目
	region_flags[0][8] = false;
	region_flags[1][8] = false;
	region_flags[2][8] = false;
	region_flags[3][8] = false;
	region_flags[4][8] = false;
	region_flags[5][8] = false;
	region_flags[6][8] = false;
	region_flags[7][8] = false;
	region_flags[8][8] = false;
	region_flags[9][8] = false;
	//10段目
	region_flags[0][9] = false;
	region_flags[1][9] = false;
	region_flags[2][9] = false;
	region_flags[3][9] = false;
	region_flags[4][9] = false;
	region_flags[5][9] = false;
	region_flags[6][9] = false;
	region_flags[7][9] = false;
	region_flags[8][9] = false;
	region_flags[9][9] = false;
	//*/
	
	int sheet_thickness = 3; // シートの厚さ（y軸方向の粒子を配置する範囲）

	for (int i = 0; i < N; ++i) {
		auto& a = sites[i];
		int unit_cell_id = site_id / 2;
		int ix = unit_cell_id % lattice_x;
		int iy = (unit_cell_id / lattice_x) % lattice_y;
		int iz = unit_cell_id / (lattice_x * lattice_y);

		if (iy < sheet_thickness) {
			// 既存の配置ロジック
			if (iz < 10) {
				// 下から1000個分の粒子はそのまま敷き詰める
				a.exist_atom_id = num_atoms;
				atoms.push_back(EventAtom{ num_atoms, site_id });
				++num_atoms;
				//敷き詰めない
				//a.exist_atom_id = SiteInfo::ATOM_NONE;
			}
			else {
				// z > 10の領域を25個のエリアに分割
				int region_x = ix / region_width_x;
				int region_z = (iz - 10) / region_height_z;

				if (region_flags[region_x][region_z] && dist(mt) < material_density) {
					a.exist_atom_id = num_atoms;
					atoms.push_back(EventAtom{ num_atoms, site_id });
					++num_atoms;
				}
				else {
					a.exist_atom_id = SiteInfo::ATOM_NONE;
				}
			}
		}
		++site_id;
	}
}

vec3d GetCoordinate(int site_id, int lattice_x, int lattice_y, int lattice_z, double lattice_constant) {
	int odd = site_id % 2;  //偶数なら0奇数なら1
	int unit_cell_id = site_id / 2;
	int ix = unit_cell_id % lattice_x;
	int iy = (unit_cell_id / lattice_x) % lattice_y;
	int iz = unit_cell_id / (lattice_x * lattice_y);
	return vec3d{ lattice_constant * ((double)ix * 2 + odd), lattice_constant * ((double)iy * 2 + odd),lattice_constant * ((double)iz * 2 + odd) };
}


//隣のサイトを洗い出す関数//
//
void GetMoveTargetList(std::vector<int>& neighbor_ids, int site_id, int lattice_x, int lattice_y, int lattice_z) {
	const int num_in_cell = 2;
	int odd = site_id % num_in_cell;  //偶数なら0奇数なら1
	int unit_cell_id = site_id / num_in_cell;
	int ix = unit_cell_id % lattice_x;
	int iy = (unit_cell_id / lattice_x) % lattice_y;
	int iz = unit_cell_id / (lattice_x * lattice_y);

	auto GetSiteID = [&](int id_in_cell, int ix, int iy, int iz) {
		ix = (ix + lattice_x) % lattice_x; //periodic//
		iy = (iy + lattice_y) % lattice_y; //periodic//
		iz = (iz + lattice_z) % lattice_z; //periodic//
		return id_in_cell + num_in_cell * (ix + lattice_x * (iy + lattice_y * iz));
		};

	if (odd == 0) {
		neighbor_ids.push_back(GetSiteID(1, ix-1, iy-1, iz));
		neighbor_ids.push_back(GetSiteID(1, ix, iy-1, iz));
		neighbor_ids.push_back(GetSiteID(1, ix-1, iy, iz));
		neighbor_ids.push_back(GetSiteID(1, ix, iy, iz));
		/*
		neighbor_ids.push_back((((ix - 1 + lattice_x) % lattice_x) + lattice_x * ((iy - 1 + lattice_y) % lattice_y)) * 2 + 1 + lattice_x * lattice_y * iz * 2);
		neighbor_ids.push_back((ix + lattice_x * ((iy - 1 + lattice_y) % lattice_y)) * 2 + 1 + lattice_x * lattice_y * iz * 2);
		neighbor_ids.push_back((((ix - 1 + lattice_x) % lattice_x) + lattice_x * iy) * 2 + 1 + lattice_x * lattice_y * iz * 2);
		neighbor_ids.push_back((ix + lattice_x * iy) * 2 + 1 + lattice_x * lattice_y * iz * 2);
		*/
		
		//neighbor_ids.push_back(GetSiteID(1, ix - 1, iy - 1, iz-1));
		//neighbor_ids.push_back(GetSiteID(1, ix, iy - 1, iz - 1));
		//neighbor_ids.push_back(GetSiteID(1, ix - 1, iy, iz - 1));
		//neighbor_ids.push_back(GetSiteID(1, ix, iy, iz - 1));
		//if (iz > 0) {
		neighbor_ids.push_back(GetSiteID(1, ix - 1, iy - 1, iz - 1));
		neighbor_ids.push_back(GetSiteID(1, ix, iy - 1, iz - 1));
		neighbor_ids.push_back(GetSiteID(1, ix - 1, iy, iz - 1));
		neighbor_ids.push_back(GetSiteID(1, ix, iy, iz - 1));

			/*
			neighbor_ids.push_back((((ix - 1 + lattice_x) % lattice_x) + lattice_x * ((iy - 1 + lattice_y) % lattice_y)) * 2 + 1 + lattice_x * lattice_y * (iz - 1) * 2);
			neighbor_ids.push_back((ix + lattice_x * ((iy - 1 + lattice_y) % lattice_y)) * 2 + 1 + lattice_x * lattice_y * (iz - 1) * 2);
			neighbor_ids.push_back((((ix - 1 + lattice_x) % lattice_x) + lattice_x * iy) * 2 + 1 + lattice_x * lattice_y * (iz - 1) * 2);
			neighbor_ids.push_back((ix + lattice_x * iy) * 2 + 1 + lattice_x * lattice_y * (iz - 1) * 2);
			*/
		//}
			
	}
	else {

		neighbor_ids.push_back(GetSiteID(0, ix + 1, iy + 1, iz));
		neighbor_ids.push_back(GetSiteID(0, ix, iy + 1, iz));
		neighbor_ids.push_back(GetSiteID(0, ix + 1, iy, iz));
		neighbor_ids.push_back(GetSiteID(0, ix, iy, iz));
		/*
		neighbor_ids.push_back((((ix + 1 ) % lattice_x) + lattice_x * ((iy + 1 ) % lattice_y)) * 2 + lattice_x * lattice_y * iz * 2);
		neighbor_ids.push_back((ix + lattice_x * ((iy + 1 ) % lattice_y)) * 2 + lattice_x * lattice_y * iz * 2);
		neighbor_ids.push_back((((ix + 1 ) % lattice_x) + lattice_x * iy) * 2 + lattice_x * lattice_y * iz * 2);
		neighbor_ids.push_back((ix + lattice_x * iy) * 2 + lattice_x * lattice_y * iz * 2);
		*/
		
		//neighbor_ids.push_back(GetSiteID(0, ix + 1, iy + 1, iz + 1));
		//neighbor_ids.push_back(GetSiteID(0, ix, iy + 1, iz + 1));
		//neighbor_ids.push_back(GetSiteID(0, ix + 1, iy, iz + 1));
		//neighbor_ids.push_back(GetSiteID(0, ix, iy, iz + 1));

		//if (iz < lattice_z - 1) {
		neighbor_ids.push_back(GetSiteID(0, ix + 1, iy + 1, iz + 1));
		neighbor_ids.push_back(GetSiteID(0, ix, iy + 1, iz + 1));
		neighbor_ids.push_back(GetSiteID(0, ix + 1, iy, iz + 1));
		neighbor_ids.push_back(GetSiteID(0, ix, iy, iz + 1));

			/*
			neighbor_ids.push_back((((ix + 1 ) % lattice_x) + lattice_x * ((iy + 1 ) % lattice_y)) * 2 + lattice_x * lattice_y * (iz + 1) * 2);
			neighbor_ids.push_back((ix + lattice_x * ((iy + 1 ) % lattice_y)) * 2 + lattice_x * lattice_y * (iz + 1) * 2);
			neighbor_ids.push_back((((ix + 1 ) % lattice_x) + lattice_x * iy) * 2 + lattice_x * lattice_y * (iz + 1) * 2);
			neighbor_ids.push_back((ix + lattice_x * iy) * 2 + lattice_x * lattice_y * (iz + 1) * 2);
			*/
		//}
		
	}

}
std::vector<int> GetNeighborAtomPresence(const std::vector<SiteInfo>& sites, int site_id, int lattice_x, int lattice_y, int lattice_z) {
	std::vector<int> neighbor_ids;
	GetMoveTargetList(neighbor_ids, site_id, lattice_x, lattice_y, lattice_z);

	std::vector<int> neighbor_presence;
	for (const auto& id : neighbor_ids) {
		if (sites[id].exist_atom_id != SiteInfo::ATOM_NONE) {
			neighbor_presence.push_back(1); // 原子が存在する場合
		}
		else {
			neighbor_presence.push_back(0); // 原子が存在しない場合
		}
	}
	return neighbor_presence;
}

int FindoutMoveTarget(EventAtom& a, const std::vector<SiteInfo>& sites, int lattice_x, int lattice_y, int lattice_z) {
	const int current_id = a.currest_site_id;

	if (sites[current_id].exist_atom_id == 2595) {
		int abc = current_id;
	}

	//粒子に設定済みの移動リストをクリア//
	a.paths.clear();

	std::vector<int> target_ids;
	GetMoveTargetList(target_ids, current_id, lattice_x, lattice_y, lattice_z);

	// 現在位置の周辺8サイトのIDを取得
	std::vector<int> current_neighbor_ids;
	GetMoveTargetList(current_neighbor_ids, current_id, lattice_x, lattice_y, lattice_z);

	//エネルギー差 E_C - E_A ////////////////////////////////////////////
	// 注目原子の周辺原子数を数えて、その平方根をE_aに加算
	int N_current = 0;
	for (const auto& tid : target_ids) {
		if (sites[tid].exist_atom_id != SiteInfo::ATOM_NONE) {
			N_current++;
		}
	}
	double dE_C_from_A = D * (double)N_current + C * sqrt((double)N_current);
	//0-(-C*rootN)

	// 周辺8サイトについて、各サイトの周辺原子数を数えて、その平方根をE_aに加算

	//点oまわりのループ
	for (const auto& neighbor_id : current_neighbor_ids) {
		if (sites[neighbor_id].exist_atom_id == SiteInfo::ATOM_NONE) continue;

		int N_neighbor = 0;
		std::vector<int> neighbor_target_ids;
		GetMoveTargetList(neighbor_target_ids, neighbor_id, lattice_x, lattice_y, lattice_z);


		for (const auto& tid : neighbor_target_ids) {
			if (sites[tid].exist_atom_id != SiteInfo::ATOM_NONE) {
				N_neighbor++;
			}
		}

		if (N_neighbor > 0) {
			dE_C_from_A +=  -D * ((double)(N_neighbor - 1) - (double)N_neighbor) - C * (sqrt((double)(N_neighbor - 1)) - sqrt((double)N_neighbor));
		}
	}

	////////////////////////////////////////////エネルギー差 E_C - E_A //

	

	int num = 0;
	for (const auto& tid : target_ids) {
		if (sites[tid].exist_atom_id == SiteInfo::ATOM_NONE) {

			//エネルギー差 E_B - E_C ////////////////////////////////////////////


			// 移動先サイトに他の粒子がいないとき
			// 行先に原子aがいると仮定したときのE_bを算出
			// 行先をbとして
			// サイトbの周りの原子を数えてN_bとする
			// E_b = -sqrt(N_1) - sqrt(N_2) - sqrt(N_3) - sqrt(N_4) - ... - sqrt(N_8) - sqrt(N_b)
			// 行き先bはtid
			// 行先bから見た周辺サイトをリストアップ
			std::vector<int> surround_b;
			GetMoveTargetList(surround_b, tid, lattice_x, lattice_y, lattice_z);



			// 注目原子の周辺原子数を数えて、その平方根をE_bに加算
			int N_target = 0;
			for (const auto& id : surround_b) {
				if (sites[id].exist_atom_id != SiteInfo::ATOM_NONE && id != current_id) {
					N_target++;
				}
			}

			double dE_B_from_C= -D * (double)N_target - C * (sqrt((double)N_target));
			

			
			
			// 周辺8サイトについて、各サイトの周辺原子数を数えて、その平方根をE_bに加算
			for (const auto& neighbor_id : surround_b) {
				if (sites[neighbor_id].exist_atom_id == SiteInfo::ATOM_NONE ) continue;
				if (neighbor_id == current_id) continue;

				int N_neighbor = 0;
				std::vector<int> neighbor_target_ids;
				GetMoveTargetList(neighbor_target_ids, neighbor_id, lattice_x, lattice_y, lattice_z);

				for (const auto& id : neighbor_target_ids) {
					if (sites[id].exist_atom_id != SiteInfo::ATOM_NONE && id != current_id) {
						N_neighbor++;
					}
				}

				//dE_B_from_C +=  -D * ((double)N_neighbor + 1 - (double)N_neighbor) - C * (sqrt((double)N_neighbor + 1) - sqrt((double)N_neighbor));
				if (N_neighbor == 0) {
					dE_B_from_C = 0;
				}
				else {
					dE_B_from_C += -D * ((double)N_neighbor + 1 - (double)N_neighbor) - C * (sqrt((double)N_neighbor + 1) - sqrt((double)N_neighbor));
				}
			}

			double dE_B_from_A = dE_B_from_C + dE_C_from_A;

			////////////////////////////////////////////エネルギー差 E_B - E_C //


			// 暫定的にkbT=1.0
			double ratio;


			if ((dE_B_from_A) > 0.0) {
				ratio = exp(-dE_B_from_A / kbT);
			}
			else {
				ratio = 1.0;
			}

			a.paths.emplace_back(MoveTarget{ current_id, tid, ratio });
			++num;
		}
	}

	return num;
}


void GetSecondNeighborList(std::vector<int>& neighbor_ids, int site_id, int lattice_x, int lattice_y, int lattice_z) {
	const int num_in_cell = 2;
	int id_in_cell = site_id % num_in_cell;
	int unit_cell_id = site_id / num_in_cell;
	int ix = unit_cell_id % lattice_x;
	int iy = (unit_cell_id / lattice_x) % lattice_y;
	int iz = unit_cell_id / (lattice_x * lattice_y);

	auto GetSiteID = [&](int id_in_cell, int ix, int iy, int iz) {
		ix = (ix + lattice_x) % lattice_x; //periodic//
		iy = (iy + lattice_y) % lattice_y; //periodic//
		iz = (iz + lattice_z) % lattice_z; //periodic//
		return id_in_cell + num_in_cell * (ix + lattice_x * (iy + lattice_y * iz));
		};


	neighbor_ids.push_back(GetSiteID(id_in_cell, ix + 1, iy, iz));
	neighbor_ids.push_back(GetSiteID(id_in_cell, ix, iy + 1, iz));
	neighbor_ids.push_back(GetSiteID(id_in_cell, ix, iy, iz + 1));
	neighbor_ids.push_back(GetSiteID(id_in_cell, ix - 1, iy, iz));
	neighbor_ids.push_back(GetSiteID(id_in_cell, ix, iy - 1, iz));
	neighbor_ids.push_back(GetSiteID(id_in_cell, ix, iy, iz - 1));

}


void GetResearchSiteList(std::vector<int>& neighbor_ids, int site_id_before, int site_id_after, int lattice_x, int lattice_y, int lattice_z) {
	

	GetMoveTargetList(neighbor_ids, site_id_before, lattice_x, lattice_y, lattice_z);
	GetMoveTargetList(neighbor_ids, site_id_after, lattice_x, lattice_y, lattice_z);
	GetSecondNeighborList(neighbor_ids, site_id_before, lattice_x, lattice_y, lattice_z);
	GetSecondNeighborList(neighbor_ids, site_id_after, lattice_x, lattice_y, lattice_z);

	std::sort(neighbor_ids.begin(), neighbor_ids.end());
	neighbor_ids.erase(std::unique(neighbor_ids.begin(), neighbor_ids.end()), neighbor_ids.end());

}


int main(int argc, char* argv[]) {
	printf("Simple KMC start--------------------\n");
	// 開始時刻を記録
	auto start_time = std::chrono::high_resolution_clock::now();
	const int64_t STEPS = 300000;
	const int lattice_x = 50;
	const int lattice_y = 2;
	const int lattice_z = 60;
	double lattice_constant = 3.0;
	double box_axis_org[12]{ lattice_constant * (double)lattice_x * 2, 0.0, 0.0,
					0.0, lattice_constant * (double)lattice_y * 2, 0.0,
					0.0, 0.0, lattice_constant * (double)lattice_z * 2,
					0.0, 0.0, 0.0 };

	std::vector<SiteInfo> sites(lattice_x * lattice_y * lattice_z * 2);
	std::vector<EventAtom> atoms;

	//データ出力機能を持ったクラス//
	KMCLogger logger;
	FILE* fp = fopen("kmc_log.krb", "w");
	logger.Initialize(fp, box_axis_org);

	//サイトに初期の粒子を配置//
	InitializeSite(sites, atoms,lattice_x, lattice_y,  lattice_z);

	//初期粒子の位置をloggerに登録//
	for (const auto& a : atoms) {
		const vec3d position = GetCoordinate(a.currest_site_id,lattice_x, lattice_y, lattice_z, lattice_constant);
		logger.Add(0, a.currest_site_id, a.currest_site_id, position, a.atom_id);
	}

	KMCEventList kmc_event_list;
	//全ての粒子の移動候補を洗い出す//
	for (auto& a : atoms) {
		const int num_found = FindoutMoveTarget(a, sites, lattice_x, lattice_y, lattice_z);
		if (num_found > 0) {
			double ratio = a.Ratio();
			const int key = a.atom_id; //粒子IDをユニークkeyとしてイベント登録する//
			kmc_event_list.Insert(key, ratio, &a);
		}
	}
	

	//乱数の準備//
	const unsigned int seed = 82109832;
	std::mt19937 mt(seed);
	std::uniform_real_distribution<double> dist(0.0,1.0);
	double elapse_time = 0.0;
	std::vector<double> time_list;
	for (int64_t istep = 1; istep <= STEPS; ++istep) {
		//(1)イベントを起こす//
		//(1.1)乱数の生成
		const double total_ratio = kmc_event_list.TotalRatio();
		const double point = dist(mt) * total_ratio;

		//(1.2)移動粒子を選択//
		double residual;
		auto hit = kmc_event_list.Bring(point, &residual);
		auto& a = hit.value();
		
		// 選択された粒子の周辺粒子数をカウント
		int neighbor_count = 0;
		std::vector<int> neighbor_ids;
		GetMoveTargetList(neighbor_ids, a->currest_site_id, lattice_x, lattice_y, lattice_z);
		for (const auto& id : neighbor_ids) {
			if (sites[id].exist_atom_id != SiteInfo::ATOM_NONE) {
				neighbor_count++;
			}
		}

		//(1.3)選択された粒子は複数の移動先を持つのでそのうち一つを選択
		double sum=0.0;
		int path_index = 0;
		for (const auto& p : a->paths) {
			if (residual <= sum + p.ratio) {
				//bring//
				break;
			} else {
				sum += p.ratio;
				++path_index;
			}
		}

		//(2)選択された粒子を移動させる//
		int site_id_before = a->currest_site_id;
		int site_id_after = a->paths[path_index].size_id_dest;
		sites[site_id_before].exist_atom_id = SiteInfo::ATOM_NONE;
		sites[site_id_after].exist_atom_id = a->atom_id;
		a->currest_site_id = site_id_after;

		//(2.2)移動をログに登録//
		const vec3d position = GetCoordinate(site_id_after, lattice_x, lattice_y, lattice_z, lattice_constant);
		logger.Add(istep, site_id_before, site_id_after, position, a->atom_id);

		//(3)移動先の再探索//
		//粒子の移動に伴い、移動した粒子の新しい移動先を見つける//
		//加えてその周辺の粒子も移動候補と移動確率が変化しているので再設定//
		//(3.1)再探索が必要な粒子をリストアップ//
		std::vector<int> re_search_atoms;
		GetResearchSiteList(re_search_atoms, site_id_before, site_id_after, lattice_x, lattice_y, lattice_z);
		for (const auto& id : re_search_atoms) {
			//idは再探索しなければいけないサイト//
			if (sites[id].exist_atom_id != SiteInfo::ATOM_NONE) {
				//再探索しなければいけないサイトに粒子がいたので再探索する//
				EventAtom& a = atoms[sites[id].exist_atom_id];
				const int num_found = FindoutMoveTarget(a, sites, lattice_x, lattice_y, lattice_z);
				if (num_found > 0) {
					double ratio = a.Ratio();
					const int key = a.atom_id; //粒子IDをユニークkeyとしてイベント登録する//
					kmc_event_list.Insert(key, ratio, &a);//updateの場合も自動で判別//
				} else {//(num_found==0)//
					//移動先が無くなったのでイベントは起こらないようにする//
					const int key = a.atom_id; //粒子IDをユニークkeyとしてイベント削除する//
					kmc_event_list.Erase(key);
				}

			}
		}
		//時間の積算
		elapse_time += 1.0 / total_ratio;
		time_list.push_back(elapse_time);
		// 10000ステップごとに情報を表示
		if (istep % 100000 == 0) {
			auto current_time = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time);

			std::cout << "Step: " << istep << ", Time: " << duration.count() << " seconds" << std::endl;

			
		}
	}
	
	//粒子の移動のログを出力
	logger.Flush(fp);
	fclose(fp);

	//経過時間の出力
	{
		FILE* fp = fopen("elapse_time.txt", "w");
		const int i_end = time_list.size();
		for (int i = 0; i < i_end; ++i) {
			fprintf(fp, "%d\t%.15f\n", i, time_list[i]);
		}
		fclose(fp);
	}
	return 0;
}
