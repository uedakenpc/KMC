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
	int size_id_org;      //�ړ��O�̃T�C�g������ID//
	int size_id_dest;     //�ړ���̃T�C�g������ID//
	double ratio;         //�ړ��p�x(�P�ʎ��ԓ�����̔����m��)//
};

//����T�C�g���猩�ė�(�ړ���)�ƂȂ�T�C�g(�p�X)�̐��̏��//
//���[�U�[�̐ݒ肷��i�q�Ɉˑ����ēK�؂Ȓl��ݒ肹��//

//constexpr size_t NUM_NEIGHBOR_SITES = 4; 
struct EventAtom {//����̗��q�ł������̈ړ�������̂ŁA��������̏��ɓZ�߂�ׂ̃N���X//
	int atom_id;            //���q�ŗL��ID
	int currest_site_id;    //���q�����݈ʒu���Ă���T�C�g//
	std::vector<MoveTarget> paths;
	//std::array<MoveTarget, NUM_NEIGHBOR_SITES> paths;
	//int num_effectivepaths = 0;    //�L���Ȉړ���̐�//

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
* ��{���[��
* - ���q�͌��߂�ꂽ�T�C�g(�i�q�_)�̏���ړ�������̂Ƃ���
* - �T�C�g�͏ꏊ���ƂɈ�ӂɌ��߂�ꂽ�ŗL��ID������
*   - �t��ID��������΃T�C�g�̍��W�Ȃǂ���������̂Ƃ���
* 
***********************************************/


struct SiteInfo {
	int exist_atom_id;  //�T�C�g�ɑ��݂��闱�q��ID//���݂��Ȃ��Ƃ���ATOM_NONE�ƂȂ�
	static const int ATOM_NONE = -1;
};

/**
* �����̗��q�z�u�����߂�
* �v�Z�������ޗ��̍\���ɍ��킹�ēK�X����������
* ���̗�ł�sites�ŗ^�����Ă���T�C�g�ɁA����m���ɏ]���ė��q��z�u����
**/
void InitializeSite(std::vector<SiteInfo>& sites, std::vector<EventAtom>& atoms) {
	unsigned int seed = 123456789;
	std::mt19937 mt(seed);
	std::uniform_real_distribution<double> dist(0.0, 1.0);

	int num_atoms = 0;
	int site_id = 0;
	for (auto&& a : sites) {
		if (dist(mt) < material_density) {//���q������//
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
	int odd = site_id % 2;  //�����Ȃ�0��Ȃ�1
	int unit_cell_id = site_id / 2;
	int ix = unit_cell_id % lattice_x;
	int iy = (unit_cell_id / lattice_x) % lattice_y;
	int iz = unit_cell_id / (lattice_x * lattice_y);
	return vec3d{ lattice_constant * (double)ix * 2 + lattice_constant * odd, lattice_constant * (double)iy * 2 + lattice_constant * odd,lattice_constant * (double)iz * 2 + lattice_constant * odd };
}

void GetMoveTargetList(std::vector<int>& neighbor_ids, int site_id, int lattice_x, int lattice_y, int lattice_z) {
	int odd = site_id % 2;  //�����Ȃ�0��Ȃ�1
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
			neighbor_presence.push_back(1); // ���q�����݂���ꍇ
		}
		else {
			neighbor_presence.push_back(0); // ���q�����݂��Ȃ��ꍇ
		}
	}
	return neighbor_presence;
}

int FindoutMoveTarget(EventAtom& a, const std::vector<SiteInfo>& sites, int lattice_x, int lattice_y, int lattice_z) {
	const int current_id = a.currest_site_id;

	//���q�ɐݒ�ς݂̈ړ����X�g���N���A//
	a.paths.clear();

	//(1)���݈ʒu����̈ړ���ƂȂ�T�C�g�����X�g�A�b�v//
	std::vector<int> target_ids;
	GetMoveTargetList(target_ids, current_id, lattice_x, lattice_y, lattice_z);
	//(2)���݈ʒu�̃G�l���M�[�v�ZE_A
	//a�̎��ӗ��q����G�l���M�[�����߂�
	//E_A=-sqrt(N_A)
	//target_id �̓��X�g�A�b�v�ς�
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
			//�ړ���T�C�g�ɑ��̗��q�����Ȃ��Ƃ�//
			//�ړ��m�������߂�//
			//�����ł͎b��I��1.0//
			//double ratio = 1.0;
			//���q�̈ړ���ƈړ��m����ǉ�//
			//ratio�̓G�l���M�[���Ō��߂�

			//�s��Ɍ��qa������Ɖ��肵���Ƃ���E_b���Z�o
			//�s���b�Ƃ���
			//�T�C�gb�̎���̌��q�𐔂���N_b�Ƃ���
			//E_b=-sqrt(N_b)
			//E_a�Ɠ�����`
			//�s����b�́@tid
			//�s��b���猩�����ӃT�C�g�����X�g�A�b�v
			std::vector<int> surround_b;
			GetMoveTargetList(surround_b, tid, lattice_x, lattice_y, lattice_z);
			//N_b�𐔂���
			int N_b = 0;
			for (const auto& id2 : surround_b) {
				if (sites[id2].exist_atom_id != SiteInfo::ATOM_NONE) {
					if (id2 != current_id) {
						N_b++;
					}
					

			

				}
			}
			//���̎��_��N_b������ꂽ
			double E_b = 0.0;
			if (N_b > 0) {
				E_b = -sqrt((double)N_b);
			}

			//�b��I��kbT=1.0
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

	//�f�[�^�o�͋@�\���������N���X//
	KMCLogger logger;
	FILE* fp = fopen("kmc_log.krb", "w");
	logger.Initialize(fp, box_axis_org);

	//�T�C�g�ɏ����̗��q��z�u//
	InitializeSite(sites, atoms);

	//�������q�̈ʒu��logger�ɓo�^//
	for (const auto& a : atoms) {
		const vec3d position = GetCoordinate(a.currest_site_id,lattice_x, lattice_y, lattice_z, lattice_constant);
		logger.Add(0, a.currest_site_id, a.currest_site_id, position, a.atom_id);
	}

	KMCEventList kmc_event_list;
	//�S�Ă̗��q�̈ړ�����􂢏o��//
	for (auto& a : atoms) {
		const int num_found = FindoutMoveTarget(a, sites, lattice_x, lattice_y, lattice_z);
		if (num_found > 0) {
			double ratio = a.Ratio();
			const int key = a.atom_id; //���qID�����j�[�Nkey�Ƃ��ăC�x���g�o�^����//
			kmc_event_list.Insert(key, ratio, &a);
		}
	}
	

	//�����̏���//
	const unsigned int seed = 82109832;
	std::mt19937 mt(seed);
	std::uniform_real_distribution<double> dist(0.0,1.0);

	for (int64_t istep = 1; istep <= STEPS; ++istep) {
		//(1)�C�x���g���N����//
		//(1.1)�����̐���
		const double total_ratio = kmc_event_list.TotalRatio();
		const double point = dist(mt) * total_ratio;

		//(1.2)�ړ����q��I��//
		double residual;
		auto hit = kmc_event_list.Bring(point, &residual);
		auto& a = hit.value();
		
		//(1.3)�I�����ꂽ���q�͕����̈ړ�������̂ł��̂������I��
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

		//(2)�I�����ꂽ���q���ړ�������//
		int site_id_before = a->currest_site_id;
		int site_id_after = a->paths[path_index].size_id_dest;
		sites[site_id_before].exist_atom_id = SiteInfo::ATOM_NONE;
		sites[site_id_after].exist_atom_id = a->atom_id;
		a->currest_site_id = site_id_after;

		//(2.2)�ړ������O�ɓo�^//
		const vec3d position = GetCoordinate(site_id_after, lattice_x, lattice_y, lattice_z, lattice_constant);
		logger.Add(istep, site_id_before, site_id_after, position, a->atom_id);

		//(3)�ړ���̍ĒT��//
		//���q�̈ړ��ɔ����A�ړ��������q�̐V�����ړ����������//
		//�����Ă��̎��ӂ̗��q���ړ����ƈړ��m�����ω����Ă���̂ōĐݒ�//
		//(3.1)�ĒT�����K�v�ȗ��q�����X�g�A�b�v//
		std::vector<int> re_search_atoms;
		GetResearchSiteList(re_search_atoms, site_id_before, site_id_after, lattice_x, lattice_y, lattice_z);
		for (const auto& id : re_search_atoms) {
			//id�͍ĒT�����Ȃ���΂����Ȃ��T�C�g//
			if (sites[id].exist_atom_id != SiteInfo::ATOM_NONE) {
				//�ĒT�����Ȃ���΂����Ȃ��T�C�g�ɗ��q�������̂ōĒT������//
				EventAtom& a = atoms[sites[id].exist_atom_id];
				const int num_found = FindoutMoveTarget(a, sites, lattice_x, lattice_y, lattice_z);
				if (num_found > 0) {
					double ratio = a.Ratio();
					const int key = a.atom_id; //���qID�����j�[�Nkey�Ƃ��ăC�x���g�o�^����//
					kmc_event_list.Insert(key, ratio, &a);//update�̏ꍇ�������Ŕ���//
				} else {//(num_found==0)//
					//�ړ��悪�����Ȃ����̂ŃC�x���g�͋N����Ȃ��悤�ɂ���//
					const int key = a.atom_id; //���qID�����j�[�Nkey�Ƃ��ăC�x���g�폜����//
					kmc_event_list.Erase(key);
				}

			}
		}
		
	}
	
	//���q�̈ړ��̃��O���o��
	logger.Flush(fp);
	fclose(fp);
	
	return 0;
}
