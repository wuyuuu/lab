#include <bits/stdc++.h>
#include <iostream>
#include <chrono>
#include <random>
#include <math.h>

using namespace std;
const int n = 30;
const int psize = 50;
double lb = -100,ub = 100;
double init_eta = 1.0;
vector<vector<double>> eta;
default_random_engine generator;
normal_distribution<double> cauchy(0, 1);
cauchy_distribution<double> gaussian(0, 1);

double func(vector<double>ele){
    double temp = 0;
    for (int i = 0; i < ele.size(); i++){
        temp += abs(ele[i]) ;
    }
    double t = 1;
    for(int i=0;i<ele.size();++i){
        t*= abs(ele[i]);
    }
    return temp + t;
}

double eval(vector<double> a) {
    double ans = 0;
    return 1.0/func(a);
}

vector<double> eval_vector(vector<vector<double>> a) {
    vector<double> ans;
    for (int i = 0; i < a.size(); ++i)
        ans.push_back(eval(a[i]));
    return ans;
}

vector<vector<double>> init_p(int psize) {
    vector<vector<double>> ans;
    while (psize--) {
        vector<double> temp;
        for (int i = 0; i < n; i++) {
            temp.push_back(lb + (rand() % 10000) / 10000.0 * (ub-lb));
        }
        ans.push_back(temp);
    }
//    eta = ans;
    return ans;
}


vector<double> mutate(vector<double> a, vector<double> eta_i) {
    for (int i = 0; i < a.size(); i++) {
        double d = eta_i[i] * (0.5 * gaussian(generator) + 0.5 * cauchy(generator));
        a[i] += d;
        a[i] = min(ub,a[i]);
        a[i] = max(lb,a[i]);
    }
    return a;
}

vector<int> robin_pos;

vector<int> robin(vector<vector<double>> p, int m = 10) {
    auto score = eval_vector(p);
    if (robin_pos.size() != p.size()) {
        robin_pos.clear();
        for (int i = 0; i < p.size(); ++i)robin_pos.push_back(i);
    };
    vector<int> ans(p.size(), 0);
    for (int i = 0; i < p.size(); i++) {
        random_shuffle(robin_pos.begin(), robin_pos.end());
        int cnt = 0;
        for (int j = 0; j < m; j++) {
            if (score[robin_pos[j]] < score[i])ans[i]++;
        }
    }
    sort(robin_pos.begin(), robin_pos.end(), [&](int x, int y) {
        return score[x] > score[y];
    });
    return robin_pos;
}

double eta_n = 1.0;

vector<double> update_eta(vector<double> e) {
    double N = gaussian(generator);
    double tau = 1.0 / sqrt(2.0 * sqrt(eta_n));
    double tau_ = 1.0 / sqrt(2.0 * eta_n);
    for (int j = 0; j < e.size(); j++) {
        e[j] *= exp(tau_ * N + tau * gaussian(generator));
        e[j] = max(e[j],0.0001);
    }
    eta_n += 1.0;
    return e;
}

vector<vector<double>> next_generation(vector<vector<double>> p) {
    int n_p = p.size();
    for (int i = 0; i < n_p; ++i) {
        p.push_back(mutate(p[i], eta[i]));
        eta.push_back(update_eta(eta[i]));
    }
    vector<int> pos = robin(p);
    vector<vector<double>> next_eta;
    vector<vector<double>> nex;
    for (int i = 0; i < n_p; i++) {
        int idx = pos[i];
        next_eta.push_back(eta[idx]);
        nex.push_back(p[idx]);
    }
    eta = next_eta;
    return nex;
}

double bestv = 1000000.0;



vector<double> get_score(vector<vector<double>> p) {
    vector<double> ans;
    for (auto ele:p) {
        ans.push_back(func(ele));
    }
    return ans;
}

double getbest(vector<vector<double>> p) {
    auto score = get_score(p);
    double maxv = 10000000.0;
    for (int i = 0; i < p.size(); i++) {
        bestv = min(bestv, score[i]);
        maxv = min(maxv, score[i]);
    }
    return maxv;
}


int main() {
    int experiment_time = 10;
    int max_g = 10000;
    vector<double> output(max_g, 0.0);
    vector<double>data;
    for (int i = 0; i < experiment_time; ++i) {
        bestv = 10000.0;
        cout << "i = " << i << endl;
        vector<vector<double>> p = init_p(psize);
        eta = vector<vector<double>>(psize,vector<double>(n,init_eta));
        for (int g = 0; g < max_g; g++) {
            p = next_generation(p);
            auto nowBest = getbest(p);
            output[g] += nowBest;
            if (g % (max_g / 10) == 0)
                cout << nowBest << endl;
        }
        data.push_back(getbest(p));
    }
    double mean = 0;
    for(auto v:data)mean+=v;mean/=data.size();
    printf("mean = %f\n",mean);
    double var = 0;
    for(auto v:data)var += (mean-v)*(mean-v);var /= data.size();
    printf("var = %f\n",var);

    return 0;
}
