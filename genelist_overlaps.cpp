#include <iostream>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

//class holding all annotations read in from gmt file in memory
struct annotations{
	
	string name, web_reference;
	vector<string> genes;
	
};

//toupper for argv[n] queries
void upstr(char *s){
	
	char  *p;
	for (p = s; *p != '\0'; p++) 
		*p = (char) toupper(*p);
		
}

//load all genesets to memory
vector<annotations> loadtable(char *filename){

	vector<annotations> pusher;
	string line;
	ifstream myfile (filename);
	if (myfile.is_open()){
		//grab each line of file
		while ( getline (myfile,line) ){
			//tokenize line, and push to string vector
			istringstream iss(line);
			string token;
			annotations organizer;
			vector<string> temp_holder;
			while(getline(iss,token,'\t')){
				if (token != "")
					temp_holder.push_back(token);
			}
			//reorganize strings and push to unordered struct
			organizer.name = temp_holder[0];
			organizer.web_reference = temp_holder[1];
			vector<string> gene_holder(temp_holder.begin()+2,temp_holder.end());
			sort(gene_holder.begin(), gene_holder.end());
			organizer.genes = gene_holder;
			pusher.push_back(organizer);
		}
		myfile.close();
	}
	return(pusher);
	
}

//create hashmap of annotation name to vector index
map<string, int> generate_mapping(vector<annotations> &anno){
	
	//for quick lookup of annotation struct index by pathway name
	map<string, int> pusher;
	int i = 0;
	for (annotations x : anno){
		pusher[x.name] = i+1;
		i++;
	}
	return(pusher);

}

//calculate the number of overlapping genes for 1 genelist against all in annotation
void calculate_overlaps(vector<annotations> &anno, int lookup_index){

	//grab geneset for comparisons against
	vector<string> lookup_genes = anno[lookup_index-1].genes;
	//iterate through annotations
	for (annotations x : anno){
		vector<string> overlapping_genes; 
		//get the intersection of genes between comparison
		set_intersection(lookup_genes.begin(), lookup_genes.end(), x.genes.begin(), x.genes.end(), back_inserter(overlapping_genes));
		if (overlapping_genes.size() > 0){
			//output descriptive data about pathways with overlap, and the number and names over overlapped genes
			cout << anno[lookup_index-1].name << "\t" << x.name << "\t" << lookup_genes.size() << "\t" << x.genes.size() << "\t" << overlapping_genes.size() << "\t";
			for (string strgene : overlapping_genes)
				cout << strgene << ",";
			cout << "\n";
		}
	}

}

//outputs all gene overlaps for N number of specified pathways against entire annotation
//usage: ./genelist_overlaps <mSigDB.gmt file OR similarily formated> <search string 1> <search string 2> ... > output.csv	
int main (int argc, char* argv[]){
	
	//load genelist to memory
	vector<annotations> genelist = loadtable(argv[1]);
	
	//generate hashmap of pathway names to genelist index;
	map<string, int> hashmap = generate_mapping(genelist);
	
	//generate report
	for (int i = 2; i < argc; i++){
		upstr(argv[i]);
		if (i == 2)
			cout << "pathway1_name\tpathways2_name\tpathway1_number\tpathway2_number\toverlap_number\toverlapped_genes\n";
		if(hashmap[argv[i]] != 0)
			calculate_overlaps(genelist, hashmap[argv[i]]);
	}
	
	return 1;
	
}
