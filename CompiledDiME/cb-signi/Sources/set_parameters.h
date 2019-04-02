

#include "./standard_package/standard_include.cpp"


void statement(char *b) {
	
	cout<<"********************************************************************"<<endl;
	//cout<<"To run the program type you need to use at least options -f and -c to specify input files"<<endl;
	cout<<"\nMANDATORY OPTIONS:"<<endl<<endl;
	cout<<"-f <filename>\n\t\tThis option is used to specify the file where the network data are stored:"<<endl;
	cout<<"\t\tfilename should include a list of edges and their weights (optionally). ";
	cout<<"Only integers number are allowed. This is the format:";
	cout<<"\n\t\tsource_node target_node (weight)\n";
	cout<<"\t\tIndeed, since the program works only on undirected graphs, it does not care about which node is the source and which one is the target. ";
	cout<<"Repetitions and self-loops are neglected"<<endl;
	cout<<"-c <filename2>\n\t\tfilename2 is the community file, i.e. it must contain the groups:"<<endl;
	cout<<"\t\teach row should contain the nodes belonging to the same community separated by a tab or a blank space."<<endl<<endl;
	cout<<"\nOTHER OPTIONS:"<<endl<<endl;
	cout<<"-t <number> :\n\t\tset the threshold equal to number"<<endl;
	cout<<"-cscore :\n\t\tthe program uses the original definition of the c-score: this is slower than the default definition which gives comparable results"<<endl;
	cout<<"-nobcore :\n\t\tthe program doesn't search for the b-q-core. this can be used if you want a fast computation of the c-cores"<<endl;
	cout<<"-list :\n\t\twhen using this option the community file is expected to be a list of \"integer integer\" which are considered node-membership"<<endl;
	cout<<"-v :\n\t\tthis flags increments the verbosity"<<endl;

	
	cout<<"\n\nExamples:\n\t"<<b<<" -f network.dat -c community.dat -t 0.01"<<endl;
	cout<<"\t"<<b<<" -f network.dat -c community_list.dat -t 0.01 -list "<<endl<<endl<<endl;


}



class Parameters {
	
	
	public:
	
		Parameters();
		~Parameters(){};
		
		string netfile;
		string comfile;
		double threshold;
		
		bool usecscore;
		bool usebqcore;
		bool uselist;
		bool verbose;
		



};

Parameters::Parameters() {
		
		
		
		
		usecscore=false;
		usebqcore=true;
		uselist=false;
		verbose=false;		
		threshold=0.05;
			
		
		
};














bool set_parameters(int argc, char * argv[], Parameters & par1) {
	
	
	

	int argct = 0;
	string temp;
	
	if (argc <= 1) { // if no arguments, return statement about program usage.
		
		statement(argv[0]);
		return false;
	}


	cout<<"*********************\tSETTING PARAMETERS\t\t*********************"<<endl;

	deque<string> command_line;
	for(int i=0; i<argc; i++) {
		string as(argv[i]);
		command_line.push_back(as);	
	}
	
	
	prints(command_line);
	int Pos=0;
	while(Pos<command_line.size()-1) {
		
		Pos++;
		
		if(command_line[Pos]=="-f") {
			
			
			Pos++;
			if(Pos==command_line.size()) {


				statement(argv[0]);
				cerr<<"ERROR reading command line..."<<endl;
				return false;
			
			}
			
			par1.netfile=command_line[Pos];
			cout<<"setting network file: "<<par1.netfile<<endl;

		
		}
		else if(command_line[Pos]=="-c") {
			
			
			Pos++;
			if(Pos==command_line.size()) {
				
				cerr<<"ERROR reading command line..."<<endl;
				return false;
			
			}
				
			par1.comfile=command_line[Pos];
			cout<<"setting community file: "<<par1.comfile<<endl;

		
		}
		else if(command_line[Pos]=="-t") {
			
			
			Pos++;

			if(Pos==command_line.size() || !cast_string_to_double(command_line[Pos], par1.threshold)) {
				
				cerr<<"ERROR reading command line..."<<endl;
				return false;
			
			}
			
			
			cout<<"setting threshold: "<<par1.threshold<<endl;

			
		
		}
		else if(command_line[Pos]=="-cscore") {
			par1.usecscore=true;
			cout<<"setting cscore option"<<endl;
		}
		
		else if(command_line[Pos]=="-nobcore") {
			par1.usebqcore=false;
			cout<<"setting no-q-bscore option"<<endl;
		}
		
		else if(command_line[Pos]=="-list") {
			par1.uselist=true;
			cout<<"setting list option"<<endl;
		}
		else if(command_line[Pos]=="-v") {
			par1.verbose=true;
			cout<<"setting high verbosity"<<endl;
		}
		else {
			
			statement(argv[0]);
			cerr<<"ERROR:: flag "<<command_line[Pos]<<" is unknown."<<endl;
			return false;
			
		}
		
		
	}
	
	
	
	if(par1.netfile.empty()) {
		
		cerr<<"network unspecified (you have to use option -f filename)"<<endl;
		return false;
	
	}
	
	if(par1.comfile.empty()) {
		
		cerr<<"community file unspecified (you have to use option -c filename)"<<endl;
		return false;
	
	}
	
	
	return true;
}




