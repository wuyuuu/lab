#include <bits/stdc++.h>
#include <iostream>
#include <chrono>
#include <random>
#include <math.h>

using namespace std;
const int n = 30;
const int psize = 30;

vector<double>uniform(vector<double>prob){
	double d = 1.0/prob.size();
	vector<double>temp;
	for(int i=0;i<prob.size();i++)temp.push_back(d);
	return temp;
} 

vector<double>rank_vector(vector<double>prob){
    vector<int>pos;
    for(int i=0;i<prob.size();++i)pos.push_back(i);
    sort(pos.begin(),pos.end(),[&](int i,int j){
        return prob[i]<prob[j];
    });
    vector<double>temp(pos.size(),0);
    for(int i=0;i<pos.size();++i){
        temp[pos[i]]=i+1.0;
    }
    return temp;
}
vector<double>geometric_ranking(vector<double>score,double alpha=0.1){
    vector<double>temp = rank_vector(score);
    double miu = score.size();
    double tot = 0;
    for(int i=0;i<temp.size();i++){
        temp[i] = alpha*pow(1-alpha,miu-1-temp[i]);
        tot += temp[i];
    }
    for(int i=0;i<temp.size();++i)temp[i]/=tot;
    return temp;
}

double eval(vector<double> a) {
    double ans = 0;
    for (int i = 0; i < a.size(); ++i) {
        ans += (i + 1) * a[i] * a[i] * a[i] * a[i];
    }
    return 1.0 / (ans + (rand() % 10000) / 10000.0);
}

double f(vector<double> a) {
    double ans = 0;
    for (int i = 0; i < a.size(); ++i) {
        ans += (i + 1.0) * a[i] * a[i] * a[i] * a[i];
    }
    return ans + (rand() % 10000) / 10000.0;
}

vector<double> denoise_eval(vector<vector<double>> a, int k) {
    vector<double> ans;
    for (int i = 0; i < a.size(); i++)
        ans.push_back(0.0);
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < a.size(); ++j)
            ans[j] += (f(a[j]));
    }
    for (int i = 0; i < a.size(); i++)ans[i] = 1.0 / ans[i];
    return ans;
}

vector<double> eval_vector(vector<vector<double>> a) {
    vector<double> ans;
    for (int i = 0; i < a.size(); ++i)
        ans.push_back(eval(a[i]));
    return ans;
}

vector<double> arith(vector<double> a, vector<double> b, double k) {
    vector<double> ans;
    for (int i = 0; i < a.size(); i++) {
        ans.push_back(a[i] * k + b[i] * (1 - k));
    }
    return ans;
}

vector<vector<double>> init_p(int psize) {
    vector<vector<double>> ans;
    while (psize--) {
        vector<double> temp;
        for (int i = 0; i < n; i++) {
            temp.push_back(-1.28 + (rand() % 10000) / 10000.0 * 2.56);
        }
        ans.push_back(temp);
    }
    return ans;
}

double gaussrand() {
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if (phase == 0) {
        do {
            double U1 = (double) rand() / RAND_MAX;
            double U2 = (double) rand() / RAND_MAX;

            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while (S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;

    return X;
}

double gaussian(double mean, double sigma) {
    return gaussrand() * sigma + mean;
}

vector<double> mutate(vector<double> a) {
    for (int i = 0; i < a.size(); i++) {
        double d = gaussian(0, 0.1);
        a[i] += d;
    }
    return a;
}

int wheel(vector<double> prob) {
    double r = rand() % 1000 * 1.0 / 1000;
    int i = 0;
    while (r > 0 and i < prob.size() - 1) {
        r -= prob[i];
        i += 1;
    }
    return i;
}

double printeval(vector<double> a) {
    double temp = 0;
    for (int i = 0; i < a.size(); i++)temp += (i + 1) * a[i] * a[i] * a[i] * a[i];
    return temp;
}

vector<vector<double>> next_generation(vector<vector<double>> p, double mutate_rate) {
    int n_p = p.size();

    vector<vector<double>> nex;

//    vector<double> score = denoise_eval(p, 33);
    vector<double> score = eval_vector(p);
    score = geometric_ranking(score);
    score = uniform(score);
    double tot = 0;
    for (auto e:score)tot += e;
    for (int i = 0; i < n_p; i++)score[i] /= tot;

    vector<int> pos;
    for (int i = 0; i < n_p; i++)pos.push_back(i);

    sort(pos.begin(), pos.end(), [&](int i, int j) {
        return score[i] > score[j];
    });

    for (int i = 0; i < n_p / 10; i++) {
        break;
        nex.push_back(p[pos[i]]);
    }

    int cnt = nex.size();
    while (cnt < 2 * n_p) {
        int i = wheel(score), j = wheel(score);
        vector<double> x = arith(p[i], p[j], 0.2), y = arith(p[j], p[i], 0.2);
        nex.push_back(x);
        nex.push_back(y);
        cnt += 2;
        if (rand() % 10000 * 1.0 / 10000.0 <= mutate_rate) {
            nex.push_back(mutate(p[i]));
            nex.push_back(mutate(p[j]));
            cnt += 2;
        }
    }
    pos.clear();
    for (int i = 0; i < cnt; i++)pos.push_back(i);
//    score = denoise_eval(nex, 33);
    score = eval_vector(nex);
//    score = geometric_ranking(score);
//    score = uniform(score);
    sort(pos.begin(), pos.end(), [&](int i, int j) {
        return score[i] > score[j];
    });
    vector<vector<double>> nexx;
    for (int i = 0; i < n_p; i++)nexx.push_back(nex[pos[i]]);
    return nexx;
}

double bestv = 1000000.0;

vector<double> get_score(vector<vector<double>> p) {
    vector<double> ans;
    for (auto ele:p) {
        double temp = 0;
        for (int i = 0; i < ele.size(); i++)
            temp += (i + 1.0) * ele[i] * ele[i] * ele[i] * ele[i];
        ans.push_back(temp);
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
    int experiment_time = 5;
    int max_g = 5000;
//    max_g /= 33;
    vector<double> output(max_g, 0.0);
    for (int i = 0; i < experiment_time; ++i) {
        bestv = 10000.0;
        cout << "i = " << i << endl;
        vector<vector<double>> p = init_p(psize);
        for (int g = 0; g < max_g; g++) {
            p = next_generation(p, 0.1);
            auto nowBest = getbest(p);
            output[g] += nowBest;
            if (g % (max_g / 10) == 0)
                cout << nowBest << endl;
        }
    }
    printf("best v = %f", bestv);
    freopen("D:\\My_Code\\EA\\Lab3_uniform.txt", "w", stdout);
    for (int i = 0; i < max_g; i++)cout << output[i] / (experiment_time * 1.0) << endl;
    return 0;
}
