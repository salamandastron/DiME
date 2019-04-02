
/*
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                               *
 *	This program is free software; you can redistribute it and/or modify         *
 *  it under the terms of the GNU General Public License as published by         *
 *  the Free Software Foundation; either version 2 of the License, or            *
 *  (at your option) any later version.                                          *
 *                                                                               *
 *  This program is distributed in the hope that it will be useful,              *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
 *  GNU General Public License for more details.                                 *
 *                                                                               *
 *  You should have received a copy of the GNU General Public License            *
 *  along with this program; if not, write to the Free Software                  *
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA    *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                               *
 *  Created by Andrea Lancichinetti on 15/5/09 (email: arg.lanci@gmail.com)      *
 *  Location: ISI foundation, Turin, Italy                                       *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */


#include "set_parameters.h"

Parameters P1;


#include "clanetwork.h"


int main(int argc, char * argv[]) {		

	srand_file();

	cout<<"erasing previous output file if any..."<<endl;
	ofstream eout1("c_core.dat");
	ofstream eout2("b_core.dat");
	ofstream eout3("c_erased.dat");
	ofstream eout4("b_erased.dat");
	
	if(set_parameters(argc, argv, P1)==false) {
		cerr<<"Read file ReadMe.txt for more information"<<endl;
		return -1;
	}

	bool good_file;
	clanetwork luca(P1.netfile, good_file);
	if(good_file==false)
		return -1;
	if(luca.size()==0) {
		
		cerr<<"network file not found"<<endl;
		return -1;
	}
	cout<<"network:: "<<luca.size()<<" nodes and "<<luca.edges()<<" edges;\t average degree = "<<2*luca.edges()/luca.size()<<endl;
	

	deque<deque<int> > ten;
	if(P1.uselist==false)
		get_partition_from_file(P1.comfile, ten);
	else
		get_partition_from_file_list(P1.comfile, ten);
	
	if(luca.translate(ten)==-1)
			return -1;
	
	
		
	cout<<ten.size()<<" cluster(s) found"<<endl;
	
	
	string sto=P1.netfile;
	sto.append(".table");
	char bsto[1000];
	cast_string_to_char(sto, bsto);
	ofstream bsoo(bsto);
	
	for(int i=0; i<ten.size(); i++) {

		deque<int> & g = ten[i];
		
		
		cout<<endl<<"================================================="<<endl;
		cout<<"Cluster "<<i+1<<":\tsize "<<g.size()<<";\tKIN: "<<luca.kin(g)<<"\tKTOT: "<<luca.ktot(g)<<endl;
		
		
		set<int> group_s;
		for (int j=0; j<g.size(); j++)
			group_s.insert(g[j]);	
		if(group_s.size()<3) {
			cerr<<"group is too small"<<endl;
		}
		else {
		
			int kw, kv_eq, degree_ref;
			int worst_node = luca.compute_the_worst(group_s, kv_eq, kw, degree_ref);
			cout<<"worst_node "<<luca.id_of(worst_node)<<";\tits internal degree "<<kw<<";\tits degree "<<degree_ref<<endl;
			
			double cscore;
			int cis;
			
			if(P1.usecscore) {
				cscore=luca.c_score(g);
				cout<<"cscore = "<<cscore<<endl;
				cis=luca.c__inside(g, P1.threshold);
				
			}
			else {
				double err;
				cscore=luca.r_score(g, err);
				cout<<"cscore "<<cscore<<" +/- "<<err<<endl;
				cis=luca.r__inside(g, P1.threshold);
			}
			
			double err;
			double bscore=luca.b_score(g, err);
			cout<<"bscore "<<bscore<<" +/- "<<err<<endl;
			if(P1.usebqcore) {
			
				int bis=luca.b__inside(g, P1.threshold);
				bsoo<<i+1<<"  "<<g.size()<<"  "<<cscore<<"  "<<bscore<<"  "<<cis<<"  "<<bis<<" "<<endl;
			
			}
			else
				bsoo<<i+1<<"  "<<g.size()<<"  "<<cscore<<"  "<<bscore<<"  "<<cis<<" "<<endl;

		}
	
	
	}
	
	
	return 0;
		
}






	
	
	
	