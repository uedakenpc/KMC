#include <cstdio>
#include <array>
#include <vector>
#include <random>
#include "AVLTreeSum3.h"
#include "KMCLogger.h"
#include <fstream>
#include <iostream>

inline static const double kbT = 0.1;

inline static const double material_density = 0.2;

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
void InitializeSite(std::vector<SiteInfo>& sites, std::vector<EventAtom>& atoms) {
	unsigned int seed = 123456789;
	std::mt19937 mt(seed);
	std::uniform_real_distribution<double> dist(0.0, 1.0);

	int num_atoms = 0;
	int site_id = 0;
	for (auto&& a : sites) {
		if (dist(mt) < material_density) {//粒子が存在//
			a.exist_atom_id = num_atoms;
			atoms.push_back(EventAtom{ num_atoms, site_id });
			++num_atoms;
		} else {
			a.exist_atom_id = SiteInfo::ATOM_NONE;
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
	return vec3d{ lattice_constant * (double)ix * 2 + lattice_constant * odd, lattice_constant * (double)iy * 2 + lattice_constant * odd,lattice_constant * (double)iz * 2 + lattice_constant * odd };
}

void GetMoveTargetList(std::vector<int>& neighbor_ids, int site_id, int lattice_x, int lattice_y, int lattice_z) {
	int odd = site_id % 2;  //偶数なら0奇数なら1
	int unit_cell_id = site_id / 2;
	int ix = unit_cell_id % lattice_x;
	int iy = (unit_cell_id / lattice_x) % lattice_y;
	int iz = unit_cell_id / (lattice_x * lattice_y);

	if (odd == 0) {
		neighbor_ids.push_back((((ix - 1 + lattice_x) % lattice_x) + lattice_x * ((iy - 1 + lattice_y) % lattice_y)) * 2 + 1 + lattice_x * lattice_y * iz * 2);
		neighbor_ids.push_back((ix + lattice_x * ((iy - 1 + lattice_y) % lattice_y)) * 2 + 1 + lattice_x * lattice_y * iz * 2);
		neighbor_ids.push_back((((ix - 1 + lattice_x) % lattice_x) + lattice_x * iy) * 2 + 1 + lattice_x * lattice_y * iz * 2);
		neighbor_ids.push_back((ix + lattice_x * iy) * 2 + 1 + lattice_x * lattice_y * iz * 2);
		if (iz > 0) {
			neighbor_ids.push_back((((ix - 1 + lattice_x) % lattice_x) + lattice_x * ((iy - 1 + lattice_y) % lattice_y)) * 2 + 1 + lattice_x * lattice_y * (iz - 1) * 2);
			neighbor_ids.push_back((ix + lattice_x * ((iy - 1 + lattice_y) % lattice_y)) * 2 + 1 + lattice_x * lattice_y * (iz - 1) * 2);
			neighbor_ids.push_back((((ix - 1 + lattice_x) % lattice_x) + lattice_x * iy) * 2 + 1 + lattice_x * lattice_y * (iz - 1) * 2);
			neighbor_ids.push_back((ix + lattice_x * iy) * 2 + 1 + lattice_x * lattice_y * (iz - 1) * 2);
		}
	}
	else {
		neighbor_ids.push_back((((ix - 1 + lattice_x) % lattice_x) + lattice_x * ((iy - 1 + lattice_y) % lattice_y)) * 2 + lattice_x * lattice_y * iz * 2);
		neighbor_ids.push_back((ix + lattice_x * ((iy - 1 + lattice_y) % lattice_y)) * 2 + lattice_x * lattice_y * iz * 2);
		neighbor_ids.push_back((((ix - 1 + lattice_x) % lattice_x) + lattice_x * iy) * 2 + lattice_x * lattice_y * iz * 2);
		neighbor_ids.push_back((ix + lattice_x * iy) * 2 + lattice_x * lattice_y * iz * 2);
		if (iz < lattice_z - 1) {
			neighbor_ids.push_back((((ix - 1 + lattice_x) % lattice_x) + lattice_x * ((iy - 1 + lattice_y) % lattice_y)) * 2 + lattice_x * lattice_y * (iz + 1) * 2);
			neighbor_ids.push_back((ix + lattice_x * ((iy - 1 + lattice_y) % lattice_y)) * 2 + lattice_x * lattice_y * (iz + 1) * 2);
			neighbor_ids.push_back((((ix - 1 + lattice_x) % lattice_x) + lattice_x * iy) * 2 + lattice_x * lattice_y * (iz + 1) * 2);
			neighbor_ids.push_back((ix + lattice_x * iy) * 2 + lattice_x * lattice_y * (iz + 1) * 2);
		}
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

	//粒子に設定済みの移動リストをクリア//
	a.paths.clear();

	//(1)現在位置からの移動先となるサイトをリストアップ//
	std::vector<int> target_ids;
	GetMoveTargetList(target_ids, current_id, lattice_x, lattice_y, lattice_z);
	//(2)現在位置のエネルギー計算E_A
	//aの周辺粒子からエネルギーを決める
	//E_A=-sqrt(N_A)
	//target_id はリストアップ済み
	int N_a = 0;
	for (const auto& tid : target_ids) {
		if (sites[tid].exist_atom_id != SiteInfo::ATOM_NONE) {
			
			N_a++;
			
		}
	}

	double E_a = 0.0;
	if (N_a > 0) {
		E_a = -sqrt((double)N_a);
	}

	
	int num=0;
	for (const auto& tid : target_ids) {
		if (sites[tid].exist_atom_id == SiteInfo::ATOM_NONE) {
			//移動先サイトに他の粒子がいないとき//
			//移動確率を決める//
			//ここでは暫定的に1.0//
			//double ratio = 1.0;
			//原子の移動先と移動確率を追加//
			//ratioはエネルギー差で決める

			//行先に原子aがいると仮定したときのE_bを算出
			//行先をbとして
			//サイトbの周りの原子を数えてN_bとする
			//E_b=-sqrt(N_b)
			//E_aと同じ定義
			//行き先bは　tid
			//行先bから見た周辺サイトをリストアップ
			std::vector<int> surround_b;
			GetMoveTargetList(surround_b, tid, lattice_x, lattice_y, lattice_z);
			//N_bを数える
			int N_b = 0;
			for (const auto& id2 : surround_b) {
				if (sites[id2].exist_atom_id != SiteInfo::ATOM_NONE) {
					if (id2 != current_id) {
						N_b++;
					}
					

			

				}
			}
			//この時点でN_bが得られた
			double E_b = 0.0;
			if (N_b > 0) {
				E_b = -sqrt((double)N_b);
			}

			//暫定的にkbT=1.0
			double ratio;
			if ((E_b - E_a) > 0.0) {
				ratio = exp(-(E_b - E_a) / kbT);
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



void GetResearchSiteList(std::vector<int>& neighbor_ids, int site_id_before, int site_id_after, int lattice_x, int lattice_y, int lattice_z) {
	

	GetMoveTargetList(neighbor_ids, site_id_before, lattice_x, lattice_y, lattice_z);
	GetMoveTargetList(neighbor_ids, site_id_after, lattice_x, lattice_y, lattice_z);

	std::sort(neighbor_ids.begin(), neighbor_ids.end());
	neighbor_ids.erase(std::unique(neighbor_ids.begin(), neighbor_ids.end()), neighbor_ids.end());

}


int main(int argc, char* argv[]) {
	printf("Simple KMC start--------------------\n");
	const int64_t STEPS = 10000;
	const int lattice_x = 10;
	const int lattice_y = 1;
	const int lattice_z = 10;
	double lattice_constant = 1.0;
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
	InitializeSite(sites, atoms);

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

	for (int64_t istep = 1; istep <= STEPS; ++istep) {
		//(1)イベントを起こす//
		//(1.1)乱数の生成
		const double total_ratio = kmc_event_list.TotalRatio();
		const double point = dist(mt) * total_ratio;

		//(1.2)移動粒子を選択//
		double residual;
		auto hit = kmc_event_list.Bring(point, &residual);
		auto& a = hit.value();
		
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
		
	}
	
	//粒子の移動のログを出力
	logger.Flush(fp);
	fclose(fp);
	
	return 0;
}
