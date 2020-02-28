//
// Created by Wuy19 on 2020/2/28.
//
#include <bits/stdc++.h>

using namespace std;
const int n = 5;

double eval(int n) {
    return 1.0 / (n * n + 1.0);
}

vector<double> eval_vector(vector<int> a) {
    vector<double>ans;
    for (int i = 0; i < a.size(); ++i)
	ans.push_back(eval(a[i]));
    return ans;
}

vector<int> init_p(int psize) {
    vector<int> ans;
    while (psize--) {
        ans.push_back((rand()%(1<<n)));
    }
    return ans;
}

int crossover(int a, int b) {
    int mask = rand() % (1 << n);
    return (a & mask) + (b & (~mask));
}

int one_bit_flipping(int a) {
    int pos = rand() % n;
    int mask = 1 << pos;
    return a ^ mask;
}

int multi_bit_flipping(int a, int k) {
    while (k--) {
        int pos = rand() % n;
        int mask = 1 << pos;
        a ^= mask;
    }
    return a;
}

int bitwsie_mutation(int a) {
    double p = 1.0 / n;
    for (int i = 0; i < n; i++) {
        if ((rand() % 10000) / 10000.0 <= p) {
            a ^= (1 << i);
        }
    }
    return a;
}

int wheel(vector<double> prob) {
    double r = rand() % 10000 * 1.0 / 10000.0;
    int i = 0;
    while (r > 0 and i < prob.size() - 1) {
        r -= prob[i];
        i += 1;
    }
    return i;
}

vector<int> next_generation(vector<int> p, double mutate_rate) {
	
    vector<double> score = eval_vector(p);
    
    double tot = 0;
    for (auto e:score)tot += e;
    for (int i = 0; i < score.size(); i++)score[i] /= tot;
    
    vector<int> nex;
    while (nex.size() < p.size()) {
        int i = wheel(score), j = wheel(score);
        int x = crossover(i, j), y = crossover(j, i);
        nex.push_back(i);nex.push_back(j);
        if (rand() % 10000 * 1.0 / 10000.0 < mutate_rate){
//            nex.push_back(one_bit_flipping(i));
//            nex.push_back(one_bit_flipping(j));
            nex.push_back(bitwsie_mutation(i));
            nex.push_back(bitwsie_mutation(j));
        }

    }
    while(nex.size() >p.size())nex.pop_back();
    return nex;
}

bool best(vector<int> p) {
    for (auto num:p)if (!num)return true;
    return false;
}

int main() {
	int experiment_time = 1000;
	int tot_g = 0;
	for(int i=0;i<experiment_time ;++i){
    		vector<int> p = init_p(5);
    		int g = 0;
    		while (!best(p)) {
        		g++;
        		p = next_generation(p,0.5);
    		}
    		tot_g += g;
   	 }
    printf("avg g = %f",tot_g*1.0/experiment_time);
    return 0;
}
