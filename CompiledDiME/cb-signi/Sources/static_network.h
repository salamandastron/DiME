#if !defined(STATIC_static_network_INCLUDED)
#define STATIC_static_network_INCLUDED


#include "./standard_package/standard_include.cpp"
#include "wsarray.h"








class static_network {
	
	
	
	public:
		
		static_network(deque<deque<int> > & , deque<int> & );
		static_network(bool , deque<deque<int> > & , deque<deque<double> > & , deque<int> & );
		static_network(string, bool &);
		~static_network();
		

		int draw(string, bool);
		int draw(string);
		int draw_consecutive(string , string , bool);
		int draw_consecutive(string , string);

		void print_id(const deque<int> & a, ostream &);
		void print_id(deque<deque<int> > & , ostream &);
		void print_id(deque<set<int> > & , ostream &);
		void print_id(set<int> & , ostream & );
		void deque_id(deque<int> & );
		
		int connected(string );		
		void print_subgraph(string , const deque<int> & );
		void set_subgraph(deque<int> & , deque<deque<int> > & , deque<deque<double> > & );
		
		void print_random_version(string );
		void print_connected_components(ostream &);
		void same_component(int , set<int> &);
		void set_connected_components(deque<deque<int> > & );
		void set_connected_components(deque<set<int> > & );
		
		int translate (deque<deque<int> > &);
		void get_id_label (map <int, int> &);
		int get_degree(deque<int> &);

		int size() {return dim;};
		double edges() {return tstrength;};
		
		double kin (const deque<int> &);
		double kin(const set<int> &);
		double ktot (const deque<int> &);
		double ktot (const set<int> &);
		
		void erase_link(int , int);
		
		
		double newman_modularity(const double kin_g, const double ktot_g);		
		double newman_modularity(set<int> & s);
		double newman_modularity(deque<set<int> > & );
		double newman_modularity(deque<int> & s);
		double newman_modularity(deque<deque<int> > & );
		
		
		int compute_betweeness_single_source(int,  map<pair<int, int>, double> &);
		int component_betweeness(int , map<pair<int, int>, double> & , set<int> & );
		int all_betweeness(map<pair<int, int>, double> & );
		
		int id_of(int a) {return vertices[a]->id_num;};
		int set_mem_adj (deque<deque<int> > & );		// unweighted networks!

		
		
	protected:




		class  vertex {
				
			public:
				
				vertex(int , int , int);
				~vertex();
						
				double kplus(const deque<int> &);
				double kplus(const set<int> &);
							
				int id_num;
				double strength;
				wsarray* links;
				

		};

				
				
		int dim;									// number of nodes
		double tstrength;							// number of links (weighted)
		bool weighted;
		
		deque <vertex*> vertices;
		
		
		
	private:
	
		void same_component__(int , set<int> & );
		void propagate_bw(int* , int* , int, set<int> & , set<int> &, deque<int> & );
			
};







static_network::vertex::vertex(int b, int c, int preall) {
	
	id_num=b;
	links=new wsarray(preall);
	
	
}



static_network::vertex::~vertex() {
	
	delete links;
	links=NULL;
	
	
	
}



double static_network::vertex::kplus(const deque<int> &a) {

	// computes the internal degree of the vertex respect with a
	
	double f=0;
	for (int i=0; i<a.size(); i++)
		f+=links->posweightof(a[i]).second;
		
		
	return f;
	
}





double static_network::vertex::kplus(const set<int> &a) {
	
	// computes the internal degree of the vertex respect with a (a is supposed to be sorted)
	
	double f=0;
	
	for (int i=0; i<links->size(); i++)
		if(a.find(links->l[i])!=a.end())
			f+=links->w[i];
		
	return f;
	
}






	
static_network::static_network(string file_name, bool & good_file) {
	
	
	
	
	
	
	char b[file_name.size()+1];
	cast_string_to_char(file_name, b);
	
	weighted=false;
	
	
	
	map<int, int> newlabels;
	deque<int> link_i;
	
	good_file=true;
	
	
	{
	
		int label=0;
		
	
		ifstream inb(b);
		
		string ins;
		
		
		while(getline(inb, ins)) {
			
			
			deque<string> ds;
			separate_strings(ins, ds);
			double inn1, inn2;
			if(ds.size()<2 || cast_string_to_double(ds[0], inn1)==false || cast_string_to_double(ds[1], inn2)==false) {
				cerr<<"From file "<<file_name<<": string not readable "<<ins<<" "<<endl;
				good_file=false;
				break;
			}
			
			
			else {
				
				
				int innum1=cast_int(inn1);
				int innum2=cast_int(inn2);
				
				
				map<int, int>::iterator itf=newlabels.find(innum1);
				if (itf==newlabels.end()) {
					newlabels.insert(make_pair(innum1, label++));
					link_i.push_back(1);
				}
				else 
					link_i[itf->second]++;
				
				
				itf=newlabels.find(innum2);
				if (itf==newlabels.end()) {
					newlabels.insert(make_pair(innum2, label++));
					link_i.push_back(1);
				}
				else
					link_i[itf->second]++;

					
				
			
			
			}
			
			
			
		}
	
	}
	
	
	dim=newlabels.size();
	
	for(int i=0; i<dim; i++)
		vertices.push_back(new vertex (0, 0, link_i[i]));
	
	
	for(map<int, int>::iterator itm=newlabels.begin(); itm!=newlabels.end(); itm++)
		vertices[itm->second]->id_num=itm->first;
	
	
	
	if(good_file) {
	
	
		ifstream inb(b);
		string ins;			
		while(getline(inb, ins)) {
			
			deque<string> ds;
			separate_strings(ins, ds);
	
		
			int innum1=cast_int(cast_string_to_double(ds[0]));
			int innum2=cast_int(cast_string_to_double(ds[1]));
					
			double w=1;
			if(ds.size()>2) {
				if(cast_string_to_double(ds[2], w)==false)
					w=1;
				else
					weighted=true;
			}
															
					
			int new1=newlabels[innum1];
			int new2=newlabels[innum2];
			
			
			if(new1!=new2) {		// no self loops!
				
				vertices[new1]->links->push_back(new2, w);
				vertices[new2]->links->push_back(new1, w);
				
			}

		}
		
		
		tstrength=0;

		for(int i=0; i<dim; i++) {
		
			vertices[i]->links->freeze();
			
			double strength_i=0;
			for(int j=0; j<vertices[i]->links->size(); j++)
				strength_i+=vertices[i]->links->w[j];
			
			vertices[i]->strength=strength_i;
			tstrength+=strength_i;
				
			
			
		}
		


		tstrength=tstrength/2.;

	}
	else
		cerr<<"File corrupted"<<endl;
	
	
	
	
	
	
}




static_network::static_network(deque<deque<int> > & link_per_node, deque<int> & label_rows) {
	
	
	
	
	// there is no check between label_rows and link per node but they need to have the same labels
	// link_per_node and weights_per_node are the list of links and weights. label_rows[i] is the label corresponding to row i
	
	
	weighted=false;
		
		
	map<int, int> newlabels;		// this maps the old labels with the new one
	for(int i=0; i<label_rows.size(); i++)
		newlabels.insert(make_pair(label_rows[i], newlabels.size()));


		
		
	dim=newlabels.size();
	
	
	
	for(int i=0; i<dim; i++)
		vertices.push_back(new vertex (0, 0, link_per_node[i].size()));
	
	
	for(map<int, int>::iterator itm=newlabels.begin(); itm!=newlabels.end(); itm++)
		vertices[itm->second]->id_num=itm->first;
	
	
		
	
			
	
	for(int i=0; i<link_per_node.size(); i++) {
			
		
		for(int j=0; j<link_per_node[i].size(); j++) {
			
			int new2=newlabels[link_per_node[i][j]];
			vertices[i]->links->push_back(new2, 1.);
			
		}
	
	}
	
	tstrength=0;

	for(int i=0; i<dim; i++) {
	
		vertices[i]->links->freeze();
		
		double strength_i=0;
		for(int j=0; j<vertices[i]->links->size(); j++)
			strength_i+=vertices[i]->links->w[j];
		
		vertices[i]->strength=strength_i;
		tstrength+=strength_i;
			
		
		
	}
	


	tstrength=tstrength/2.;

}



static_network::static_network(bool weighted_, deque<deque<int> > & link_per_node, deque<deque<double> > & weights_per_node, deque<int> & label_rows) {
	
	
	// there is no check between label_rows and link per node but they need to have the same labels
	// link_per_node and weights_per_node are the list of links and weights. label_rows[i] is the label corresponding to row i
	
	weighted=weighted_;
		
		
	map<int, int> newlabels;		// this maps the old labels with the new one
	for(int i=0; i<label_rows.size(); i++)
		newlabels.insert(make_pair(label_rows[i], newlabels.size()));


		
		
	dim=newlabels.size();
	
	
	
	for(int i=0; i<dim; i++)
		vertices.push_back(new vertex (0, 0, link_per_node[i].size()));
	
	
	for(map<int, int>::iterator itm=newlabels.begin(); itm!=newlabels.end(); itm++)
		vertices[itm->second]->id_num=itm->first;
	
	
		
	
			
	
	for(int i=0; i<link_per_node.size(); i++) {
			
		
		for(int j=0; j<link_per_node[i].size(); j++) {
			
			int new2=newlabels[link_per_node[i][j]];
			vertices[i]->links->push_back(new2, weights_per_node[i][j]);
			
		}
	
	}
	
	tstrength=0;

	for(int i=0; i<dim; i++) {
	
		vertices[i]->links->freeze();
		
		double strength_i=0;
		for(int j=0; j<vertices[i]->links->size(); j++)
			strength_i+=vertices[i]->links->w[j];
		
		vertices[i]->strength=strength_i;
		tstrength+=strength_i;
			
		
		
	}
	





	tstrength=tstrength/2.;

}








int static_network::get_degree(deque<int> &d) {

	d.clear();
	for (int i=0; i<dim; i++)
		d.push_back(vertices[i]->links->size());



	return 0;
}





static_network::~static_network() {
	
	for (int i=0; i<vertices.size(); i++) {
		
		delete vertices[i];
		vertices[i]=NULL;
	
	}


}




double static_network::newman_modularity(const double kin_g, const double ktot_g) {
	
	return ((kin_g)/(2.*tstrength)-pow((ktot_g)/(2.*tstrength),2));

}
		
double static_network::newman_modularity(set<int> & s) {

	return newman_modularity(kin(s), ktot(s));

};



double static_network::newman_modularity(deque<set<int> > & Comps) {

	double mm=0;
	for(int i=0; i<Comps.size(); i++)
		mm+=newman_modularity(Comps[i]);
	
	return mm;

}


double static_network::newman_modularity(deque<int> & s) {

	return newman_modularity(kin(s), ktot(s));

};



double static_network::newman_modularity(deque<deque<int> > & Comps) {

	double mm=0;
	for(int i=0; i<Comps.size(); i++)
		mm+=newman_modularity(Comps[i]);
	
	return mm;

}




void static_network::erase_link(int a, int b) {


	double sa =vertices[a]->links->posweightof(b).second;	
	
	vertices[a]->links->erase(b);
	vertices[b]->links->erase(a);

	vertices[a]->strength-=sa;
	vertices[b]->strength-=sa;

	
	tstrength-=sa;



}





int static_network::connected (string str) {

	int spannet=0;
	deque <set<int> > partic;
	deque <int> present;
	present.assign(dim,0);

	
	while (spannet!=dim) {
		
		set <int> connected;
		set <int> newcon;
				

		for (int i=0; i<dim; i++)
			if (present[i]==0) {
				connected.insert(i);
				newcon.insert(i);
				present[i]=1;
				break;
			}
			
		
	
		while (newcon.size()!=0) {
			
			set <int> nnewcon=newcon;
			newcon.clear();
			set <int>::iterator it=nnewcon.begin();
			while (it!=nnewcon.end()) {
				
				int near=0; 
				
				while (near!=vertices[*it]->links->size()) {
					present[*it]=1;
					if (connected.insert(vertices[*it]->links->l[near]).second)
						newcon.insert(vertices[*it]->links->l[near]);
					near++;
				}
				it++;
			}
		}
		
		partic.push_back(connected);
		spannet+=connected.size();
	}
	
	
	char B[100];
	cast_string_to_char(str, B);
	ofstream con_out(B);

	//cout<<"number of connected components = "<<partic.size()<<endl;
	//cout<<"dimensions"<<endl;
	
	
	int max_pcon=0;

	for (int i=0; i<partic.size(); i++) {
		
		//cout<<partic[i].size()<<"   ";
		
		
		if (partic[i].size()>=partic[max_pcon].size())
			max_pcon=i;
	}
	
	
	//cout<<endl<<endl;
	
	
	
	
	set <int>::iterator it=partic[max_pcon].begin();
	
		
	
	it=partic[max_pcon].begin();
	while (it!=partic[max_pcon].end()) {
		
		for (int j=0; j<vertices[*it]->links->size(); j++) if(vertices[*it]->id_num < vertices[vertices[*it]->links->l[j]]->id_num)
				con_out<<vertices[*it]->id_num<<"\t"<<vertices[vertices[*it]->links->l[j]]->id_num<<"\t"<<vertices[*it]->links->w[j]<<endl;
		
		it++;
	}

	return (partic.size());
}


double static_network::kin (const deque<int> & seq) {
	
	double k=0;
	for (int i=0; i<seq.size(); i++)
		k+=vertices[seq[i]]->kplus(seq);

	return k;

}


double static_network::ktot (const deque<int> &seq) {
	
	double k=0;
	for (int i=0; i<seq.size(); i++)
		k+=vertices[seq[i]]->strength;
	return k;

}


double static_network::ktot(const set <int> &s) {
	
	double k=0;
	for (set<int>::iterator it=s.begin(); it!=s.end(); it++)
		k+=vertices[*it]->strength;
	
	
	
	return k;

}


double static_network::kin(const set <int> &s) {
	
	double k=0;
	for (set<int>::iterator it=s.begin(); it!=s.end(); it++)
		k+=vertices[*it]->kplus(s);
	
	
	
	return k;

}







int static_network::draw(string file_name) {

	return draw(file_name, weighted);

}

int static_network::draw_consecutive(string file_name1, string file_name2) {
	
	return draw_consecutive(file_name1, file_name2, weighted);
}

int static_network::draw_consecutive(string file_name1, string file_name2, bool _weighted_) {
	
	
	
	
	char b[128];
	cast_string_to_char(file_name1, b);
	
	
	cout<<"drawing in file "<<b<<endl;
	ofstream graph_out(b);
	
	
	
	if (_weighted_) {
		
		for (int i=0; i<vertices.size(); i++)
			for (int j=0; j<vertices[i]->links->size(); j++) if(vertices[i]->id_num < vertices[vertices[i]->links->l[j]]->id_num)
				graph_out<<i<<"\t"<<vertices[i]->links->l[j]<<"\t"<<vertices[i]->links->w[j]<<endl;
	
	}
	
	else {
		
		for (int i=0; i<vertices.size(); i++)
			for (int j=0; j<vertices[i]->links->size(); j++) if(vertices[i]->id_num < vertices[vertices[i]->links->l[j]]->id_num)
				graph_out<<i<<"\t"<<vertices[i]->links->l[j]<<endl;
	
	
	}
	
	char bb[128];

	cast_string_to_char(file_name2, bb);
	ofstream graph_out2(bb);
	for (int i=0; i<vertices.size(); i++)
		graph_out2<<i<<" "<<vertices[i]->id_num<<endl;

		

	return 0;

}

int static_network::draw(string file_name, bool _weighted_) {
	
	
	
		int h= file_name.size();
		
		char b[h+1];
		for (int i=0; i<h; i++)
			b[i]=file_name[i];
		b[h]='\0';
		
		
		
		ofstream graph_out(b);

		
		if (_weighted_) {
			
			for (int i=0; i<vertices.size(); i++)
				for (int j=0; j<vertices[i]->links->size(); j++) if(vertices[i]->id_num < vertices[vertices[i]->links->l[j]]->id_num)
					graph_out<<vertices[i]->id_num<<"\t"<<vertices[vertices[i]->links->l[j]]->id_num<<"\t"<<vertices[i]->links->w[j]<<endl;
		
		}
		
		else {
			
			for (int i=0; i<vertices.size(); i++)
				for (int j=0; j<vertices[i]->links->size(); j++) if(vertices[i]->id_num < vertices[vertices[i]->links->l[j]]->id_num)
					graph_out<<vertices[i]->id_num<<"\t"<<vertices[vertices[i]->links->l[j]]->id_num<<endl;
		
		
		}

	return 0;

}


void static_network::print_random_version(string file_name) {
	
	
	char b[100];
	cast_string_to_char(file_name, b);
	ofstream outt(b);
	
	deque<int> nodes;
	deque<int> degrees;
	
	for(int i=0; i<dim; i++)
		nodes.push_back(vertices[i]->id_num);
	
	for(int i=0; i<dim; i++)
		degrees.push_back(cast_int(vertices[i]->strength));


	/*
	cout<<"nodes"<<endl;
	prints(nodes);
	
	cout<<"degrees"<<endl;
	prints(degrees);
	//*/
	

	
		
	// this function is to build a network with the labels stored in nodes and the degree seq in degrees (correspondence is based on the vectorial index)
	// the only complication is that you don't want the npdes to have neighbors they already have
	
	
	
	// labels will be placed in the end
	deque<set<int> > en; // this is the E of the subgraph

	{
		set<int> first;
		for(int i=0; i<nodes.size(); i++) 
			en.push_back(first);
	}
	
	
	
	multimap <int, int> degree_node;
	
	for(int i=0; i<degrees.size(); i++)
		degree_node.insert(degree_node.end(), make_pair(degrees[i], i));
	
	int var=0;

	while (degree_node.size() > 0) {
		
		multimap<int, int>::iterator itlast= degree_node.end();
		itlast--;
		
		multimap <int, int>::iterator itit= itlast;
		deque <multimap<int, int>::iterator> erasenda;
		
		int inserted=0;
		
		for (int i=0; i<itlast->first; i++) {
			
			if(itit!=degree_node.begin()) {
			
				itit--;
				
				
				en[itlast->second].insert(itit->second);
				en[itit->second].insert(itlast->second);
				inserted++;
				
				erasenda.push_back(itit);				
				
			}
			
			else
				break;
		
		}
		
		
		for (int i=0; i<erasenda.size(); i++) {
			
			
			if(erasenda[i]->first>1)
				degree_node.insert(make_pair(erasenda[i]->first - 1, erasenda[i]->second));
	
			degree_node.erase(erasenda[i]);
		
		}

		
		var+= itlast->first - inserted;
		degree_node.erase(itlast);
		
	}

	
	
	// this is to randomize the subgraph -------------------------------------------------------------------
	
	for(int node_a=0; node_a<degrees.size(); node_a++) for(int krm=0; krm<en[node_a].size(); krm++) {
	
					
				
		int random_mate=irand(degrees.size()-1);
		while (random_mate==node_a)
			random_mate=irand(degrees.size()-1);
				
		
		if (en[node_a].insert(random_mate).second) {
			
			deque <int> out_nodes;
			for (set<int>::iterator it_est=en[node_a].begin(); it_est!=en[node_a].end(); it_est++) if ((*it_est)!=random_mate)
				out_nodes.push_back(*it_est);
						
										
					
			int old_node=out_nodes[irand(out_nodes.size()-1)];
					
			en[node_a].erase(old_node);
			en[random_mate].insert(node_a);
			en[old_node].erase(node_a);

										
			deque <int> not_common;
			for (set<int>::iterator it_est=en[random_mate].begin(); it_est!=en[random_mate].end(); it_est++)
				if ((old_node!=(*it_est)) && (en[old_node].find(*it_est)==en[old_node].end()))
					not_common.push_back(*it_est);
					
						
			int node_h=not_common[irand(not_common.size()-1)];
			
			en[random_mate].erase(node_h);
			en[node_h].erase(random_mate);
			en[node_h].insert(old_node);
			en[old_node].insert(node_h);
			
			
		}
		
		
	}
	
	
	for(int i=0; i<en.size(); i++)
		for(set<int>::iterator it= en[i].begin(); it!=en[i].end(); it++)
			outt<<vertices[i]->id_num<<" "<<vertices[*it]->id_num<<endl;
		

	
	

}



void static_network::get_id_label (map <int, int> &a) {
	
	for (int i=0; i<dim; i++)
		a.insert(make_pair(vertices[i]->id_num, i));


}

void static_network::deque_id(deque<int> & a) {
	
	for(int i=0; i<a.size(); i++)	
		a[i]=vertices[a[i]]->id_num;


}

void static_network::print_id(const deque<int> & a, ostream & pout) {
	
	for (int i=0; i<a.size(); i++)
		pout<<vertices[a[i]]->id_num<<"\t";
	pout<<endl;


}


void static_network::print_id(set<int> & a, ostream & pout) {
	
	for (set<int>::iterator its=a.begin(); its!=a.end(); its++)
		pout<<vertices[*its]->id_num<<"\t";
	pout<<endl;


}



void static_network::print_id(deque<deque<int> > & a, ostream & pout) {
	
	for(int i=0; i<a.size(); i++)
		print_id(a[i], pout);
	


}

void static_network::print_id(deque<set<int> > & a, ostream & pout) {
	
	for(int i=0; i<a.size(); i++)
		print_id(a[i], pout);
	


}



int static_network::translate(deque<deque<int> > & ten) {

	map<int, int> A;
	get_id_label(A);
	
	for(int i=0; i<ten.size(); i++) {
		for(int j=0; j<ten[i].size(); j++) {
			
			map<int, int>::iterator itf=A.find(ten[i][j]);
			if(itf==A.end()) {
				
				cerr<<"ERROR: the nodes in the communities are different from those ones in the network!"<<endl;
				return -1;
			
			
			}
			
			ten[i][j]=itf->second;
			
		}
		
	
	}
	
	return 0;


}




void static_network::print_subgraph(string file_name, const deque<int> & group) {


	int h= file_name.size();
	char b[h+1];
	cast_string_to_char(file_name, b);
	ofstream subout(b);

	
	
	for(int i=0; i<group.size(); i++) {
		
		int nodei=group[i];
		
		for(int j=0; j<group.size(); j++) {
			
			
			double wij=vertices[nodei]->links->posweightof(group[j]).second;
			if(wij>0 && vertices[nodei]->id_num < vertices[group[j]]->id_num)
				subout<<vertices[nodei]->id_num<<" "<<vertices[group[j]]->id_num<<" "<<wij<<endl;
			
		
		
		}
		
	
	
	
	}


}



void static_network::set_subgraph(deque<int> & group, deque<deque<int> > & link_per_node, deque<deque<double> > & weights_per_node) {
	
	
	// in this function I'm not using id's... because I want to work with the same labels (don't want to translate)
	
	
	
	
	
	sort(group.begin(), group.end());

	weights_per_node.clear();
	link_per_node.clear();
	
	for(int i=0; i<group.size(); i++) {
		
		int nodei=group[i];
		
		
		deque<int> link_i;
		deque<double> weight_i;
		
		
		for(int j=0; j<vertices[nodei]->links->size(); j++) if(binary_search(group.begin(), group.end(), vertices[nodei]->links->l[j])) {
			
			link_i.push_back(vertices[nodei]->links->l[j]);
			weight_i.push_back(vertices[nodei]->links->w[j]);
			
		}
		
		
		link_per_node.push_back(link_i);
		weights_per_node.push_back(weight_i);
		
		
	}


}





void static_network::same_component__(int source, set<int> & mates) {

	
	deque<int> new_vertices;
	
	
	for(int i=0; i<vertices[source]->links->size(); i++) {
		
		int neigh=vertices[source]->links->l[i];
		
		
		if(mates.find(neigh)==mates.end()) {		// new vertex!!!
			
			
			mates.insert(neigh);		
			new_vertices.push_back(neigh);
		
		
		} 		
	
	
	}

	

	for(int i=0; i<new_vertices.size(); i++)
		same_component__(new_vertices[i], mates);
	
		

}



void static_network::same_component(int source, set<int> & mates) {

	
	
	mates.clear();
	mates.insert(source);
	same_component__(source, mates);
	
		

}

void static_network::set_connected_components(deque<deque<int> > & comps) {

	
	
	comps.clear();
	set<int> not_assigned;
	for(int i=0; i<dim; i++)
		not_assigned.insert(i);
	
	while(not_assigned.size()>0) {
	
		
		int source = *not_assigned.begin();
		
		set<int> mates;
		same_component(source, mates);
		
		
		deque<int> ccc;
		for(set<int>::iterator its=mates.begin(); its!=mates.end(); its++) {
			ccc.push_back(*its);
			not_assigned.erase(*its);
			
		}
		
		comps.push_back(ccc);
		
		
	
	
	
	
	}
		



}



void static_network::set_connected_components(deque<set<int> > & comps) {

	
	
	comps.clear();
	set<int> not_assigned;
	for(int i=0; i<dim; i++)
		not_assigned.insert(i);
	
	while(not_assigned.size()>0) {
	
		
		int source = *not_assigned.begin();
		
		set<int> mates;
		same_component(source, mates);
		
		for(set<int>::iterator its=mates.begin(); its!=mates.end(); its++)
			not_assigned.erase(*its);
			
		
		comps.push_back(mates);
		
	
	}
	

}

void static_network::print_connected_components(ostream & outb) {

	
	
	
	int spannet=0;
	deque <set<int> > partic;
	deque <int> present;
	present.assign(dim,0);

	
	while (spannet!=dim) {
		
		set <int> connected;
		set <int> newcon;
				

		for (int i=0; i<dim; i++)
			if (present[i]==0) {
				connected.insert(i);
				newcon.insert(i);
				present[i]=1;
				break;
			}
			
		
	
		while (newcon.size()!=0) {
			
			set <int> nnewcon=newcon;
			newcon.clear();
			set <int>::iterator it=nnewcon.begin();
			while (it!=nnewcon.end()) {
				
				int near=0; 
				
				while (near!=vertices[*it]->links->size()) {
					present[*it]=1;
					if (connected.insert(vertices[*it]->links->l[near]).second)
						newcon.insert(vertices[*it]->links->l[near]);
					near++;
				}
				it++;
			}
		}
		
		partic.push_back(connected);
		spannet+=connected.size();
	}
	
	
	for(int i=0; i<partic.size(); i++) {
		for(set<int>::iterator its=partic[i].begin(); its!=partic[i].end(); its++)
			outb<<vertices[*its]->id_num<<" ";
		
		outb<<endl;
	}
	



}





void static_network::propagate_bw(int* distances, int* weights, int source, set<int> & mates, set<int> & not_leaves, deque<int> & next_shell) {

	
	mates.insert(source);
	
	for(int i=0; i<vertices[source]->links->size(); i++) {
		
		int neigh=vertices[source]->links->l[i];
		
		
		
		if(mates.find(neigh)==mates.end()) {		// new vertex
			
			
			//cout<<"new "<<vertices[neigh]->id_num<<endl;
			distances[neigh]=distances[source]+1;
			weights[neigh]=weights[source];
			mates.insert(neigh);
			next_shell.push_back(neigh);
			not_leaves.insert(source);			
		
		
		} else if(distances[neigh]==distances[source]+1) {
			
			weights[neigh]+=weights[source];
			not_leaves.insert(source);

		}
		
	
	
	}


}






int static_network::compute_betweeness_single_source(int source, map<pair<int, int>, double> & tot_edge_bw) {
	
	
	
	// this function compute the edge betweenness of all the vertices in the component of source (only for source)
	set<int> mates;
	set<int> not_leaves;
	
	
	
	
	
	int distances[dim];
	int weights[dim];
	
	
	
	
	distances[source]=0;
	weights[source]=1;
	deque<int> present_shell;
	present_shell.push_back(source);
	
	while(true) {
		
		
		if(present_shell.empty())
			break;
		
		deque<int> next_shell;
		
		for(int i=0; i<present_shell.size(); i++)
			propagate_bw(distances, weights, present_shell[i], mates, not_leaves, next_shell);
		
		//cout<<"------------ next "<<endl;
		//print_id(next_shell, cout);

		present_shell=next_shell;
		
	}
	
	/*
	prints(distances, dim);
	prints(weights, dim);
	*/
	
	
	
	//print_id(mates, cout);
	
	deque<int> leaves;
	for(set<int>::iterator its=mates.begin(); its!=mates.end(); its++)
		if(not_leaves.find(*its)==not_leaves.end())
			leaves.push_back(*its);
	
	
	//cout<<"leaves"<<endl;
	//print_id(leaves, cout);
	
	map<pair<int, int>, double> edge_bw;	// map edge-betweenness
	for(set<int>::iterator its=mates.begin(); its!=mates.end(); its++) {
		
		for(int i=0; i<vertices[*its]->links->size(); i++) if(*its < vertices[*its]->links->l[i])
			edge_bw.insert(make_pair(make_pair(*its, vertices[*its]->links->l[i]), 0));
	
	
	}
	
	
	multimap<int, int> distance_not_leaves;
	for(set<int>::iterator its=not_leaves.begin(); its!=not_leaves.end(); its++)
		distance_not_leaves.insert(make_pair(-distances[*its], *its));
	
	
	for(deque<int>::iterator its=leaves.begin(); its!=leaves.end(); its++) {
		
		for(int i=0; i<vertices[*its]->links->size(); i++) if(distances[*its]>distances[vertices[*its]->links->l[i]]) {
			
			pair<int, int> ed;
			ed.first=min(*its, vertices[*its]->links->l[i]);
			ed.second=max(*its, vertices[*its]->links->l[i]);
			
			
			
			edge_bw[ed]=double(weights[vertices[*its]->links->l[i]])/weights[*its];
			tot_edge_bw[ed]+=double(weights[vertices[*its]->links->l[i]])/weights[*its];
			
			
		
		}
			
		
	
	}
	
	
	
	for(multimap<int, int>::iterator itm=distance_not_leaves.begin(); itm!=distance_not_leaves.end(); itm++) {
		
		
		
		//cout<<"node:::: "<<itm->second<<endl;
		
		
		double sum_of_weight=0;
		
		for(int i=0; i<vertices[itm->second]->links->size(); i++) {
			
			pair<int, int> ed;
			ed.first=min(itm->second, vertices[itm->second]->links->l[i]);
			ed.second=max(itm->second, vertices[itm->second]->links->l[i]);
			
			sum_of_weight+= edge_bw[ed];
		
		}
		
		
		for(int i=0; i<vertices[itm->second]->links->size(); i++) if(distances[itm->second]>distances[vertices[itm->second]->links->l[i]]) {
			
			pair<int, int> ed;
			ed.first=min(itm->second, vertices[itm->second]->links->l[i]);
			ed.second=max(itm->second, vertices[itm->second]->links->l[i]);
			
			edge_bw[ed]=double(weights[vertices[itm->second]->links->l[i]])/weights[itm->second]*(1 + sum_of_weight);
			tot_edge_bw[ed]+=double(weights[vertices[itm->second]->links->l[i]])/weights[itm->second]*(1 + sum_of_weight);
			
			//cout<<"pred--> "<<vertices[itm->second]->links->l[i]<<" "<<double(weights[vertices[itm->second]->links->l[i]])/weights[itm->second]*(1 + sum_of_weight)<<" "<<sum_of_weight<<endl;
			
			
		}
		
		
	
	
	}
	
	
	//cout<<"************************"<<endl;
	
	
		
	
	return 0;
	
		
}



int static_network::component_betweeness(int source, map<pair<int, int>, double> & tot_edge_bw, set<int> & mates) {		// this compute the betweenness of the edges in the component of source
	
	
	
	mates.clear();
	
	same_component(source, mates);
	
	
	for(set<int>::iterator its=mates.begin(); its!=mates.end(); its++) {
		
		for(int j=0; j<vertices[*its]->links->size(); j++) if(*its < vertices[*its]->links->l[j]) {
			
			
			pair<int, int> ed;
			ed.first=*its;
			ed.second=vertices[*its]->links->l[j];
			
			tot_edge_bw[ed]=0;			
			
			
		}
	
	}



	for(set<int>::iterator its=mates.begin(); its!=mates.end(); its++)		// this must be made n times
		compute_betweeness_single_source(*its, tot_edge_bw);								// this requires m operations
	
	
	
	
	//for(map<pair<int, int>, double>::iterator itm=tot_edge_bw.begin(); itm!=tot_edge_bw.end(); itm++)
		//cout<<vertices[itm->first.first]->id_num<<" "<<vertices[itm->first.second]->id_num<<" "<<itm->second<<endl;
	
	
	
	
	return 0;


}





int static_network::all_betweeness(map<pair<int, int>, double> & tot_edge_bw) {

	// this compute the betweenness of all the edges
	
	
	
	
	set<int> not_assigned;
	for(int i=0; i<dim; i++)
		not_assigned.insert(i);
	
	while(not_assigned.size()>0) {
	
		
		int source = *not_assigned.begin();
		set<int> mates;
		
		component_betweeness(source, tot_edge_bw, mates);		
		
		for(set<int>::iterator its=mates.begin(); its!=mates.end(); its++)
			not_assigned.erase(*its);

	
	}
	
	
	//for(map<pair<int, int>, double>::iterator itm=tot_edge_bw.begin(); itm!=tot_edge_bw.end(); itm++)
		//cout<<vertices[itm->first.first]->id_num<<" "<<vertices[itm->first.second]->id_num<<" "<<itm->second<<endl;
	

		
	return 0;


}





int static_network::set_mem_adj (deque<deque<int> > & mem_adj) {		// unweighted networks!
	
	mem_adj.clear();
	
	for (int i=0; i<dim; i++) {
		
		deque <int> first(dim);
		for (int j=0; j<dim; j++)
			first[j]=0;
		
		for (int j=0; j<vertices[i]->links->size(); j++)
			first[vertices[i]->links->l[j]]=1;
		
		mem_adj.push_back(first);
		
	
	}
	
	return 0;
}










#endif




