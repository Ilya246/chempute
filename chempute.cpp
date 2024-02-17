#include <cmath>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <limits>
#include <numeric>
#include <regex>
#include <set>
#include <unordered_set>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

void splitString(const string& str, vector<string>& vec, char delim) {
    size_t last = 0;
    for (size_t i = 0; i < str.size(); i++) {
        if (str[i] == delim) {
            if (i > last) {
                vec.push_back(str.substr(last, i - last));
            }
            last = i + 1;
        }
    }
    if (last < str.size()) {
        vec.push_back(str.substr(last, str.size() - last));
    } else if (str.size() == 0) {
        vec.push_back("");
    }
}

size_t countLeadingSpaces(const string_view& str) {
    size_t scount = 0;
    for (const char& c : str) {
        if (c == ' ') {
            ++scount;
        } else {
            break;
        }
    }
    return scount;
}

void stripSpaces(string& str) {
    size_t offset = 0;
    for (size_t i = 0; i < str.size(); ++i) {
        if (offset > 0) {
            str[i - offset] = str[i];
        }
        if (str[i] == ' ') {
            ++offset;
        }
    }
    if (offset > 0) {
        str.erase(str.size() - offset);
    }
}
string cleanString(string& str) {
    size_t offset = 0;
    for (size_t i = 0; i < str.size(); ++i) {
        if (offset > 0) {
            str[i - offset] = str[i];
        }
        if (str[i] < 32) {
            ++offset;
        }
    }
    if (offset > 0) {
        str.erase(str.size() - offset);
    }
    return str;
}

struct reagent;

struct recipe_input {
    double amount;
    bool catalyst;
};

struct recipe {
    string name;
    unordered_map<string, recipe_input*> inputs;
    unordered_map<string, double> outputs;
};

struct reagent {
    string name;
    vector<recipe*> recipes_input;
    vector<recipe*> recipes_output;
};

enum parse_state {p_reactants, p_products, p_none};

unordered_map<string, reagent> reagent_map;
unordered_map<string, recipe> reaction_map;
unordered_map<string, size_t> scale_cache;

unordered_set<string> ignore_reagents;

const short ind_num = 3;
const char ind_char = ' ';
bool doScale = true;

size_t getScale(const reagent& reag, const unordered_set<string>& traversed_in) {
    if (scale_cache.contains(reag.name)) {
        return scale_cache[reag.name];
    }
    if (reag.recipes_output.empty() || !doScale || ignore_reagents.contains(reag.name)) {
        scale_cache[reag.name] = 1;
        return 1;
    }
    unordered_set<string> traversed(traversed_in);
    traversed.insert(reag.name);
    size_t min_scale = numeric_limits<size_t>::max();
    for (recipe* reci_ptr : reag.recipes_output) {
        recipe reci = *reci_ptr;
        vector<double> outs;
        for (const auto& [in_reag_name, in_param_ptr] : reci.inputs) {
            recipe_input in_param = *in_param_ptr;
            reagent in_reag = reagent_map[in_reag_name];
            if (in_reag.recipes_output.empty() || traversed.contains(in_reag.name)) {
                continue;
            }
            outs.push_back(lcm((long long)in_param.amount, (long long)getScale(in_reag, traversed)) / in_param.amount);
        }
        size_t scale = 1;
        for (double d : outs) {
            scale = lcm((long long)scale, (long long)round(d));
        }
        scale *= reci.outputs[reag.name];
        min_scale = min(scale, min_scale);
    }
    scale_cache[reag.name] = min_scale;
    return min_scale;
}
size_t getScale(const reagent& reag) {
    return getScale(reag, unordered_set<string>());
}
void analyzeReagent(const reagent& reag, unordered_map<string, double>* amounts_in, short indent_n, double req_out, const unordered_set<string>& traversed_in) {
    unordered_set<string> traversed(traversed_in);
    traversed.insert(reag.name);
    string indent(indent_n * ind_num, ind_char);
    for (recipe* reci_ptr : reag.recipes_output) {
        unordered_map<string, double>* amounts = (amounts_in == nullptr ? new unordered_map<string, double> : amounts_in);
        recipe reci = *reci_ptr;
        double scale = getScale(reag, traversed);
        if (req_out > 0) {
            scale *= req_out / scale;
        }
        scale /= reci.outputs[reag.name];
        for (const auto& [out_reag_name, out_amount] : reci.outputs) {
            reagent out_reag = reagent_map[out_reag_name];
            cout << indent << "[O] " << out_reag.name << ": " << out_amount * scale << endl;
        }
        for (const auto& [in_reag_name, in_param_ptr] : reci.inputs) {
            recipe_input in_param = *in_param_ptr;
            reagent in_reag = reagent_map[in_reag_name];
            double required_amount = in_param.catalyst ? getScale(in_reag) * ceil(in_param.amount / (double)getScale(in_reag)) : in_param.amount * scale;
            cout << indent << (in_param.catalyst ? "[C] " : "[I] ") << in_reag.name << ": " << required_amount << endl;
            if (in_reag.recipes_output.empty() || traversed.contains(in_reag.name) || ignore_reagents.contains(in_reag.name)) {
                (*amounts)[in_reag_name] = (*amounts)[in_reag_name] + required_amount;
                continue;
            }
            analyzeReagent(in_reag, amounts, indent_n + 1, required_amount, traversed);
        }
        if (amounts_in == nullptr) {
            for (const auto& [reag_name, amount] : *amounts) {
                cout << "[B] " << reag_name << ": " << amount << endl;
            }
            delete amounts;
        }
    }
}
void analyzeReagent(const reagent& reag, short indent_n, double req_out) {
    analyzeReagent(reag, nullptr, indent_n, req_out, unordered_set<string>());
}

int main(int argc, char* argv[]) {
    string ignore_path = "chempute_ignore.txt";
    string path_begin = "Reactions/";
    bool analyzeAll = false;
    if (argc > 1) {
        for (int i = 0; i < argc; i++) {
            string arg(argv[i]);
            if (arg[0] != '-' || arg.length() < 2) {
                continue;
            }
            /* if (arg[1] == '-') {
                if (arg.rfind("--parg1", 0) == 0) {
                    targetRadius = std::stod(arg.substr(8));
                } else if (arg.rfind("--parg2", 0) == 0) {
                    tickCap = std::stoi(arg.substr(7));
                } else {
                    cout << "Unrecognized argument '" << arg << "'." << endl;
                }
                continue;
            } */
            switch (arg[1]) {
                case 'p': {
                    path_begin = arg.substr(2);
                    break;
                }
                case 's': {
                    doScale = !doScale;
                    break;
                }
                case 'a': {
                    analyzeAll = !analyzeAll;
                    break;
                }
                default: {
                    cout << "Unrecognized argument '" << arg << "'." << endl;
                    break;
                }
            }
        }
    }
    ifstream ignore_stream(ignore_path);
    string to_ignore;
    while (getline(ignore_stream, to_ignore)) {
        cleanString(to_ignore);
        ignore_reagents.insert(to_ignore);
    }
    for (const auto& file : filesystem::recursive_directory_iterator(path_begin)) {
        recipe* cur_recipe = nullptr;
        reagent* cur_reagent = nullptr;
        parse_state state = p_none;
        ifstream in_stream(file.path().relative_path());
        string line;
        while (getline(in_stream, line)) {
            vector<string> split;
            splitString(line, split, ':');
            int leading_spaces = countLeadingSpaces(split[0]);
            for (string& s : split) {
                stripSpaces(s);
            }
            switch (leading_spaces) {
                case 2: {
                    state = p_none;
                    if (split[0] == "id") {
                        recipe reac;
                        reac.name = split[1];
                        reaction_map[split[1]] = reac;
                        cur_recipe = &reaction_map[split[1]];
                    } else if (split[0] == "reactants") {
                        state = p_reactants;
                    } else if (split[0] == "products") {
                        state = p_products;
                    }
                    break;
                }
                case 4: {
                    if (state == p_none || !cur_recipe) {
                        break;
                    }
                    if (!reagent_map.contains(split[0])) {
                        reagent_map[split[0]] = reagent{split[0], {}, {}};
                        if (!analyzeAll) {
                            cout << "Parsed: " << split[0] << endl;
                        }
                    }
                    cur_reagent = &reagent_map[split[0]];
                    if (state == p_products) {
                        (*cur_reagent).recipes_output.push_back(cur_recipe);
                        (*cur_recipe).outputs[cur_reagent->name] = stod(split[1]);
                    }
                    break;
                }
                case 5:
                case 6: {
                    if (state != p_reactants || !cur_reagent || !cur_recipe) {
                        break;
                    }
                    if (split[0] == "catalyst") {
                        (*cur_recipe).inputs[cur_reagent->name]->catalyst = cleanString(split[1]) == "true";
                        break;
                    }
                    if (split[0] != "amount") {
                        break;
                    }
                    (*cur_reagent).recipes_input.push_back(cur_recipe);
                    (*cur_recipe).inputs[cur_reagent->name] = new recipe_input{stod(split[1]), false};
                    break;
                }
                default: {
                    break;
                }
            }
        }
    }
    if (analyzeAll) {
        for (const auto& [name, reag] : reagent_map) {
            cout << "\n[R] " << reag.name << endl;
            analyzeReagent(reag, 0, -1);
        }
        return 0;
    }
    cout << "Parsed reactions and reagents" << endl;
    while (true) {
        cout << "Reagent to find recipe for: " << endl;
        string analyze_ID;
        cin >> analyze_ID;
        reagent reag = reagent_map[analyze_ID];
        size_t scale = getScale(reag);
        cout << "Batch size: " << scale << endl;
        cout << "Batches to attempt to make: " << endl;
        double batches = 1;
        cin >> batches;
        analyzeReagent(reag, 0, scale * batches);
    }
}
