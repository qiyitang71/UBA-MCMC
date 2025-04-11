// Copyright (c) 2025-  The ltl2pba Authors
//
// This file is a part of ltl2pba
//
// ltl2pba is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ltl2pba is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
#include <spot/parseaut/public.hh>
#include <spot/twaalgos/isunamb.hh>
#include <spot/twaalgos/hoa.hh>
#include <spot/twaalgos/sccfilter.hh>
#include <spot/twaalgos/determinize.hh>
#include <spot/twaalgos/translate.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twa/bddprint.hh>
#include <spot/misc/version.hh>
#include <spot/misc/optionmap.hh>
#include <spot/twaalgos/split.hh>
#include <spot/tl/parse.hh>
#include <spot/tl/print.hh>
#include <spot/tl/simplify.hh>


#include <bddx.h>

#include <fstream>
#include <iostream>
// #include <boost/process.hpp>

#include <ctime>
#include <map>

#include "dcamin.hh"

// namespace bp = boost::process;

void print_usage(std::ostream &os)
{
    os << "Usage: uba2pba [OPTION...] [FILENAME...]\n";
}

void print_help()
{
    print_usage(std::cout);
    std::cout <<
        R"(The tool reduces Unambiguous Büchi auomtata (UBA) into equivalent probabilistic Büchi automata (PBA).

By default, it reads a UBA from standard input
and converts it into PBA.

Input options:
    -a FILENAME reads the UBA from FILENAME instead of stdin
Output options:
    -o FILENAME writes the output from FILENAME instead of stdout

Pre- and Post-processing:
    --preprocess[=0|1]       simplify the input formula
    --merge-transitions      merge transitions in the output automaton

Miscellaneous options:
  -h, --help    print this help
  --version     print program version
)";
}

void check_cout()
{
    std::cout.flush();
    if (!std::cout)
    {
        std::cerr << "ltl2pba: error writing to standard output\n";
        exit(2);
    }
}

void output_file(spot::twa_graph_ptr aut, const char *file)
{
    const char *opts = nullptr;
    std::ofstream outfile;
    std::string file_name(file);
    outfile.open(file_name);

    spot::print_hoa(outfile, aut, opts);
    outfile.close();
}

// bool execute(std::vector<std::string> &args)
// {
//     bp::ipstream out;
//     bp::child c(args, bp::std_out > out);

//     // now we parse the output

//     for (std::string line; c.running() && std::getline(out, line);)
//     {
//         std::cout << "Output: " << line << std::endl;
//     }
//     c.wait();
//     return false;
// }

// ================= translate to unambiguous Buchi automata ================
spot::twa_graph_ptr
read_from_file(spot::bdd_dict_ptr dict, const std::string &file_name)
{

    spot::automaton_stream_parser parser(file_name);
    spot::parsed_aut_ptr parsed_aut = parser.parse(dict);
    if (parsed_aut->format_errors(std::cerr))
    {
        std::runtime_error("File " + file_name + " is not in valid HOA format");
        return nullptr;
    }
    return parsed_aut->aut;
}

// spot::twa_graph_ptr
// translate_to_uba(spot::bdd_dict_ptr dict, spot::formula f)
// {
//     std::cout << "Translating formula " << f << " to UBA" << std::endl;

//     spot::tl_simplifier simp;
//     f = simp.simplify(f);
//     spot::twa_graph_ptr aut;
//     const std::string file_name = "uba.hoa";

//     bool ltl2tgba = false;
//     if (ltl2tgba) {
//         std::vector<std::string> cmd;
//         cmd.push_back("/usr/local/bin/ltl2tgba");
//         cmd.push_back("-U");
//         cmd.push_back("-b");
//         cmd.push_back("-f");
//         cmd.push_back(spot::str_psl(f, true));
//         cmd.push_back("-o");
//         cmd.push_back(file_name);

//         if (execute(cmd))
//         {
//             std::cerr << "Error happened" << std::endl;
//         }

//         auto aut = read_from_file(dict, file_name);
//     }else {
//         spot::translator trans;
//         trans.set_pref(spot::postprocessor::Unambiguous);
//         trans.set_type(spot::postprocessor::Buchi);
//         trans.set_level(spot::postprocessor::optimization_level::Medium);
//         aut = trans.run(&f);
//         output_file(aut, file_name.c_str());
//     }
//     std::cout << "UBA has been output to " << file_name << std::endl;
//     return aut;
// }

int compute_max_outdegree(spot::twa_graph_ptr aut)
{
    int max_degree = 1;
    for (int state = 0; state < aut->num_states(); state++)
    {
        std::map<int, int> degree_map;
        for (auto &edge : aut->out(state))
        {
            int bdd_id = edge.cond.id();
            auto it = degree_map.find(bdd_id);
            if (it != degree_map.end())
            {
                degree_map[bdd_id]++;
            }
            else
            {
                degree_map[bdd_id] = 1;
            }
        }
        for (const auto &pair : degree_map)
        {
            max_degree = std::max(max_degree, pair.second);
        }
    }
    return max_degree;
}
// precondition: number > 1
int get_num_bits(int number)
{
    if (number <= 0)
        return 0;
    if (number <= 1 || number <= 2)
        return 1;
    number -= 1; // from 0 to number - 1
    int num_bits = 0;
    while (number != 0)
    {
        number >>= 1;
        num_bits++;
    }
    return num_bits;
}

// we assume that only one proposition holds
std::string
compute_bit_repr(bdd cond, const std::map<std::string, bdd>& origs)
{
    for (const auto& item: origs) {
        bdd inter = cond & item.second;
        if (inter != bdd_false())
        {   
            return item.first;
        }
    }
    throw std::runtime_error("Unknown proposition");
}

bdd
compute_bit_repr(std::string ap_name, const std::map<std::string, bdd>& new_props)
{
    bdd prop = bdd_true();
    for (const auto& item: new_props) {
        if (item.first == ap_name)
        {   
            prop &= item.second;
        }else {
            prop &= bdd_not(item.second);
        }
    }
    return prop;
}

bdd compute_bit_repr(int number, std::vector<bdd> &aux_prop, int num_bits)
{
    if (num_bits <= 0)
        return bdd_true();
    // number of bits
    int i = 0;
    int head = 1;
    bdd prop = bdd_true();
    while (i < num_bits)
    {
        if (head & number)
        {
            prop &= aux_prop[i];
        }
        else
        {
            prop &= bdd_not(aux_prop[i]);
        }
        i++;
        head <<= 1;
    }
    return prop;
}

struct CustomComparator {
    bool operator()(const std::pair<int, bool>& a, const std::pair<int, bool>& b) const {
        // Define your custom order here.
        // Example: Sort primarily by `bool` (false before true), then by `int`.
        if (a.first != b.first) {
            return a.first < b.first; // false comes before true
        }
        return a.second < b.second; // if bools are equal, compare ints
    }
};
// determinise with auxiliary propositions
spot::twa_graph_ptr
determinise_with_aux_props(spot::twa_graph_ptr aut, int max_degree)
{
    // int num_aux_props = get_num_bits(max_degree);
    // std::cout << "Num of auxiliary bits: " << num_aux_props << std::endl;

    auto dict = spot::make_bdd_dict();
    // We need to use the same dictionary
    spot::twa_graph_ptr result = make_twa_graph(dict);
    // result->copy_ap_of(aut);
    std::map<std::string, bdd> orig_props;
    std::map<std::string, bdd> new_props;
    std::set<std::string> new_ap_names;
    for (auto& item : aut->get_dict()->var_map) {
        // std::cout << item.first << " -> " << item.second << std::endl;
        orig_props.emplace(item.first.ap_name(), bdd_ithvar(item.second));
        for (int i = 0; i < max_degree; i++)
        {
            // all propositions have unique names
            std::string new_ap_name = item.first.ap_name() + "_" + std::to_string(i);
            bdd prop = bdd_ithvar(result->register_ap(new_ap_name));
            new_props.emplace(new_ap_name, prop);
            new_ap_names.insert(new_ap_name);
        }
    }
    
    // for (auto& item : result->get_dict()->var_map) {
    //     std::cout << "f=" << item.first << " var=" << item.second << std::endl;
    // }
    result->new_states(aut->num_states() + 1);
    result->set_acceptance(1, spot::acc_cond::fin({0}));
    unsigned sink_acc_state = aut->num_states();
    for (unsigned state = 0; state < aut->num_states(); state++)
    {
        std::map<int, std::set<std::pair<int, bool>>> succ_map;
        for (auto &edge : aut->out(state))
        {
            int bdd_id = edge.cond.id();
            auto it = succ_map.find(bdd_id);
            bool acc = edge.acc.has(0);
            if (it != succ_map.end())
            {
                // found the key
                it->second.insert(std::pair<int, bool>(edge.dst, acc));
            }
            else
            {
                std::set<std::pair<int, bool>> succs;
                succs.insert(std::pair<int, bool>(edge.dst, acc));
                succ_map.emplace(bdd_id, succs);
            }
        }
        // now we resolve nondeterminism
        std::set<std::string> present_ap_names;
        for (const auto &set_pair : succ_map)
        {
            bdd cond = bdd_from_int(set_pair.first);
            int number = 0;            
            std::set<std::pair<int, bool>, CustomComparator> decrease_set;
            for (const auto &pair : set_pair.second)
            {
                decrease_set.insert(pair);
            }
            for (const auto &pair : decrease_set)
            {
                std::string ap = compute_bit_repr(cond, orig_props);
                std::string new_ap = ap + "_" + std::to_string(number);
                bdd aux_cond = compute_bit_repr(new_ap, new_props);
                // seen as co-Buchi
                if (pair.second)
                {
                    result->new_acc_edge(state, pair.first, aux_cond);
                }
                else
                {
                    result->new_edge(state, pair.first, aux_cond);
                }
                number++;
                present_ap_names.insert(new_ap);
            }
        }
        // now we redirect all left conditions to sink_acc_state
        std::set<std::string> remaining_ap_names;
        for (const auto& ap_name : new_ap_names) {
            if (present_ap_names.find(ap_name) == present_ap_names.end()) {
                remaining_ap_names.insert(ap_name);
            }
        }
        for (const auto& ap_name : remaining_ap_names) {
            bdd aux_cond = compute_bit_repr(ap_name, new_props);
            result->new_edge(state, sink_acc_state, aux_cond);
        }
    }
    result->set_init_state(aut->get_init_state_number());
    result->new_edge(sink_acc_state, sink_acc_state, bdd_true());
    // for ()
    return result;
}

// ================= translate to unambiguous Buchi automata ================

void output_hoa(spot::twa_graph_ptr aut, std::string file_name)
{
    std::ofstream outfile;
    outfile.open(file_name);
    spot::print_hoa(outfile, aut);
    outfile.close();
}

// we assume that all prosition holds exclusively
int main(int argc, char *argv[])
{
    // Declaration for input options. The rest is in cola.hpp
    // as they need to be included in other files.
    std::vector<spot::twa_graph_ptr> ubas;
    bool  check = false;

    auto dict = spot::make_bdd_dict();
    spot::option_map om;

    bool use_simulation = false;
    bool merge_transitions = false;
    bool debug = false;
    bool aut_type = false;

    std::string output_filename = "";
    std::string owl_path = "";

    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "--simplify-input")
            om.set("preprocess", true);
        else if (arg == "--merge-transitions")
            merge_transitions = true;
        else if (arg == "--level")
        {
        }
        else if (arg == "-o")
        {
            if (argc < i + 1)
            {
                std::cerr << "uba2pba: Option -o requires an argument.\n";
                return 1;
            }
            else
            {
                std::string str(argv[i + 1]);
                output_filename = str;
                i++;
            }
        }else if (arg == "-a")
        {
            if (argc < i + 1)
            {
                std::cerr << "uba2pba: Option -o requires an argument.\n";
                return 1;
            }
            else
            {
                std::string str(argv[i + 1]);
                auto aut = read_from_file(dict, str);
                ubas.push_back(aut);
                i++;
            }
        }
        else if ((arg == "--help") || (arg == "-h"))
        {
            print_help();
            check_cout();
            return 0;
        }else if (arg == "-c") {
            check = true;
        }
    }

    if (ubas.empty())
    {
        std::cerr << "uba2pba: No UBA to process? "
                     "Run 'uba2pba --help' for more help.\n";
        print_usage(std::cerr);
        return 1;
    }

    for (const auto& uba : ubas)
    {
        uba->purge_dead_states();
        int max_degree; 
        // every letter has only one proposition that holds
        spot::twa_graph_ptr aut = spot::split_edges(uba);
        if (aut->get_acceptance() != spot::acc_cond::acc_code::buchi())
        {
            std::cerr << "uba2pba: needs to input a UBA\n";
            return 1;
        }

        if (spot::is_deterministic(aut))
        {
            std::cout << "The UBA is deterministic" << std::endl;
            // we need to make this automaton as co-Buchi
            max_degree = 0;
            aut = spot::dualize(aut);
        }
        else if (spot::is_unambiguous(aut))
        {
            std::cout << "The UBA is unambiguous but not deterministic" << std::endl;
            // now we need to determinise this automaton into
            // deterministic one via adding extra propositions
            // 1. we first split the edges
            // 2. compute the maximal out degree
            max_degree = compute_max_outdegree(aut);
             // 3. now determinise with extra propositions
            aut = determinise_with_aux_props(aut, max_degree);
        }else {
            throw std::runtime_error("Needs a UBA as input");
        }
        std::cout << "The max outgoing degree is " << max_degree << std::endl;
        // output_hoa(aut, "uba.hoa");
       
        std::cout << "Minimising DCA" << std::endl;
        // merge those edges
        aut->merge_edges();
        // now we call the minimisation tool
        aut->purge_dead_states();
        std::cout << "UBA has " << aut->num_states() << " states" << std::endl;

        // output_hoa(aut, "dca.hoa");

        dca_minimiser dcamin(aut);
        auto res = dcamin.minimize();
        std::cout << "Now output DCAs" << std::endl;
        std::cout << "GfG-NCA has " << res->num_states() << " states" << std::endl;
        // aut = read_from_file(dict, "mingfgca.hoa");
        if (check ) {
            std::cout << "Equivalence: " << is_equivalent(aut, res) << std::endl;
        }
        if (!merge_transitions) aut = spot::split_edges(res);
        if (output_filename != "")
        {
            output_hoa(aut, output_filename);
        }
        else
        {
            spot::print_hoa(std::cout, aut);
            std::cout << std::endl;
        }
    }

    check_cout();

    return 0;
}
