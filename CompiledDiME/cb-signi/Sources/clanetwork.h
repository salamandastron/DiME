

#include "static_network.h"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include "border_score.h"


int convert_random_variables(deque<double> & cum1, deque<double> & cum2, int nn) {

	// convert nn from distribution1 in a realization from distribution2
	
	
	
	double a1=0;
	if(nn>0)
		a1=cum1[nn-1];
	
	double a2=cum1[nn];
	
	double r= ran4() * (a2 -a1) + a1;
	r=min(cum2[cum2.size()-1], r);
	r=max(0., r);
	
	/*
	cout<<"-------------------------------"<<endl;
	prints(cum1);
	prints(cum2);
	cout<<nn<<" "<<(lower_bound(cum2.begin(), cum2.end(), r) - cum2.begin())<<" "<<r<<endl;
	*/
	return (lower_bound(cum2.begin(), cum2.end(), r) - cum2.begin());
	

}


double compute_r_hyper(int kin_node, int kout_g, int tm, int degree_node) {

	double r=1;
	if(kin_node>0)
		r=gsl_cdf_hypergeometric_Q(kin_node-1, kout_g, tm - kout_g, degree_node);
	
	
	//cout<<"RRR "<<r<<endl;
	// kin_node -1 take account of the equal in the >=
		
	return r;

}


class clanetwork : public static_network {
	
	
	public:
	
		
		
		clanetwork(string a, bool & b) : static_network(a, b){};
		~clanetwork(){};
		
		int random_group (const int , deque<int> &);
		int good_group (const int , deque<int> &, double);
		int good_group_h (const int , deque<int> &, double);
		int good_group_greedy(const int, deque<int> &);
	
		
		double modularity(double  , double );
		double modularity(deque<int> &);

		double c_score(const deque<int> &);
		double c_score(const deque<int> & group_ggg, int & worst_node);
		int compute_mcout(const set<int> & , int , int , int , deque<deque<double> > & , map<int, int> &);
		int compute_mcout_averaged(const set<int> & group_s, int degree_ref);

		int compute_the_worst(const set<int> & group, int & kv_eq, int & kw, int & degree_ref);

		double random_cscore(const set<int> & group_s);
		void compute_rs(const set<int> & , deque<pair<double, int> > & );
		double r_score(const deque<int> &, double &);
		double r_score(const set<int> &, double &);
		double r_score(const deque<int> &);
		double r_score(const set<int> &);
		double compute_score_border(set<int> & group_s, deque<int> & border_nodes);
		double comp_border(set<int> , double &);
		void compute_rs_external(const set<int> & group, deque<int> & border_nodes, deque<pair<double, int> > & rs);
		double b_score(const deque<int> &, double &);
		double b_score(const set<int> &, double &);
		double b_score(const set<int> &);
		int c__inside(const deque<int> & group_ggg, double threshold);
		int just_the_worst(const set<int> & g);
		double c_score(set<int> group_s);
		int r__inside(const deque<int> & group_d, double threshold);
		int b__inside(const deque<int> & group_d, double threshold);
		
};






int clanetwork::compute_the_worst(const set<int> & group, int & kv_eq, int & kw, int & degree_ref) {
	
	
	
	double kin_g=kin(group);
	double ktot_g=ktot(group);
	int kout_g= cast_int(ktot_g - kin_g);
	int tm = cast_int(2 * tstrength - ktot_g);
	
	
	multimap<double, int> r_node;

	for (set<int>::iterator its=group.begin(); its!=group.end(); its++) {
		
		
		int kin_i=cast_int(vertices[*its]->kplus(group));
		int strength_i=cast_int(vertices[*its]->strength);
		
		
		//cout<<"------> "<<kin_i<<" "<<kout_g + 2 * kin_i - strength_i<<" "<<tm + strength_i<<" "<<strength_i<<endl;
		double r= - compute_r_hyper(kin_i, kout_g + 2 * kin_i - strength_i, tm + strength_i, strength_i);
		
		r_node.insert(make_pair(r, *its));
	}

	
	int worst_node = r_node.begin()->second;
	
	multimap<double, int> ::iterator ittm=r_node.begin();
	ittm++;
	kw= cast_int(vertices[worst_node]->kplus(group));
	degree_ref=cast_int(vertices[worst_node]->strength);
	
	double r_threshold = - ittm->first;
	{
		
		
		int kin_i=cast_int(vertices[ittm->second]->kplus(group));
		int strength_i=cast_int(vertices[ittm->second]->strength);
		
		r_threshold =compute_r_hyper(kin_i, kout_g + 2 * kw - degree_ref, tm + strength_i, strength_i);
		
	}
	
	
	//cout<<"threshold ... "<<r_threshold<<" "<< - ittm->first<<endl;

	r_threshold=min(- ittm->first, r_threshold);
	
	

	deque<double> cum_ref_opposite;
	kout_g += 2 * kw - degree_ref;
	tm+=degree_ref;
	
	for (int i=0; i<=degree_ref; i++)
		cum_ref_opposite.push_back(gsl_cdf_hypergeometric_Q(i, kout_g, tm - kout_g, degree_ref));

	
	kv_eq=0;
	while(cum_ref_opposite[kv_eq]/r_threshold>=1)
		kv_eq++;
	
	

	return worst_node;

}


int clanetwork::compute_mcout(const set<int> & group_s, int degree_ref, int kout_g, int tm, deque<deque<double> > & cums, map<int, int> & degree_cumlabel) {
	
	
		
	
	
	
	// you can do this if you are sure that degree_ref is included in degree_hist (you are because the worst node was erased)
	deque<double> & cum_ref=cums[degree_cumlabel[degree_ref]];
		
	
	
	int mcout=0;
	for(int i=0; i<dim; i++) if(group_s.find(i)==group_s.end())
		mcout+=convert_random_variables(cums[degree_cumlabel[cast_int(vertices[i]->strength)]], cum_ref, cast_int(vertices[i]->kplus(group_s)));
	
	
	
	//cout<<"kout - mcout "<<kout_g<<" "<<mcout<<endl;	
	//cout<<"cum_ref"<<endl;
	//prints(cum_ref);
	
	
	return mcout;

}



int clanetwork::compute_mcout_averaged(const set<int> & group_s, int degree_ref) {
 
	double kin_g = kin(group_s);
	double ktot_g = ktot(group_s);
	int kout_g = cast_int(ktot_g - kin_g);
	int tm = cast_int(2 * tstrength - ktot_g);

	
	
	
	set<int> degree_hist;
	for(int i=0; i<dim; i++) if(group_s.find(i)==group_s.end())
		degree_hist.insert(cast_int(vertices[i]->strength));
	
	
	
	
	deque<deque<double> > cums;
	map<int, int> degree_cumlabel;
	for(set<int>::iterator itm=degree_hist.begin(); itm!=degree_hist.end(); itm++) {
		
		deque<double> cum2;
		degree_cumlabel.insert(make_pair(*itm, cums.size()));
			
		for (int i=0; i<=*itm; i++)
			cum2.push_back(gsl_cdf_hypergeometric_P(i, kout_g, tm - kout_g, *itm));
		
		cums.push_back(cum2);
	
	}

	
	
	
	
	
	//------------------------------------
	
	int mcout;
	deque<double> mcs;

	double avedevmax=0.1;
	
	
	int counter=0;
	
	
	while(true) {
		
		counter++;
		int mcout=compute_mcout(group_s, degree_ref, kout_g, tm, cums, degree_cumlabel);
		mcs.push_back(mcout);
		double avedev;
		if(mcs.size()>10)
			avedev=sqrt(variance_func(mcs)/mcs.size());
		else
			avedev=1;
		
		//cout<<counter<<" iteration "<<avedev<<endl;
		//prints(mcs);
		if(avedev<avedevmax)
			break;
		
		if(counter>500)
			break;
		
	
	}
	
	
	//cout<<"average "<<(cast_int(average_func(mcs)))<<endl;
	return (cast_int(average_func(mcs)));


}


double clanetwork::c_score(set<int> group_s) {
	
	
	if(group_s.size()<3) {
		return 1;
	}
	
	
		
	int kw, kv_eq, degree_ref;
	int worst_node = compute_the_worst(group_s, kv_eq, kw, degree_ref);
	
		
	if(kv_eq<kw) {
		
		cerr<<"ERROR; please report this message to andrea.lancichinetti@isi.it, thank you!"<<endl;
		return -1;
	
	}
	
	
	group_s.erase(worst_node);
	// ******************* from now on, the worst node does not belong to the cluster ***********************************
	
		
	
	int mcout= compute_mcout_averaged(group_s, degree_ref);
	
	
		

	double cscore;
	
	
	if(kw==0)
		cscore=1;
	else if (kv_eq * (dim - group_s.size()) < mcout)
		return 1;
	else {
				
		double prob=gsl_cdf_hypergeometric_P(kw-1, mcout, kv_eq * (dim - group_s.size()) - mcout, kv_eq);
		cscore=1 - pow(prob, (dim - group_s.size()));
	
	}
	
	if(cscore<0)
		cscore=0;
	
	
	
		
	
	
	return cscore;


}



double clanetwork::c_score(const deque<int> & group_ggg) {
	
	
	set<int> group_s;
	
	for (int i=0; i<group_ggg.size(); i++)
		group_s.insert(group_ggg[i]);

		
	
	return c_score(group_s);


}









int clanetwork::random_group(const int size, deque<int> & group) {



	group.clear();
	
	if(size>dim) {
	
		cerr<<"size of the group bigger than dim"<<endl;
		return -1;
		
	}
	
	
	
	
	{
		set<int> g;
		
		while (g.size()<size)
			g.insert(irand(dim-1));
		
		
		
		for(set<int>::iterator it=g.begin(); it!=g.end(); it++)
			group.push_back(*it);
	}
	

	return 0;



}




double clanetwork::modularity(double  kin_g, double  ktot_g) {
	
	
	
	return ((kin_g)/(2.*tstrength)-pow((ktot_g)/(2.*tstrength),2));


}


double clanetwork::modularity(deque<int> &group) {
	
	
	
	double kin_g=kin(group);
	double ktot_g=ktot(group);
	

	return modularity(kin_g, ktot_g);

}


// in this version self-links are not expected
int clanetwork::good_group(const int size, deque<int> & group, double initial_temp) {
	
	
	group.clear();
	
	if (size==dim) {
	
		cerr<<"please, it is all the network!"<<endl;
		return -1;
		
	}
	
	random_group(size, group);
	
	
	
	double kin_g=kin(group);
	double ktot_g=ktot(group);
	
	double qtot= modularity(kin_g, ktot_g);
	
	double kin_i[dim];
	for (int i=0; i<dim; i++)
		kin_i[i]=vertices[i]->kplus(group);
	
	set<int> group_s;
	for (int i=0; i<group.size(); i++)
		group_s.insert(group[i]);
		
	
		
	double tstep=0.9999;
	int max_loops=3000;
	
	
	double temp=initial_temp;
	double qprevious=qtot;
	int stopper=0;
	
	int counter=0;

	while (true) {
		
		temp=temp*tstep;
		qprevious=qtot;
		
		for (int h=0; h<dim; h++) {
			
			counter++;
			int new_candidate=irand(dim-1);
			while(group_s.find(new_candidate)!=group_s.end()) 
				new_candidate=irand(dim-1);
			
			
			int old_position=irand(group.size()-1);
			int old_candidate=group[old_position];
				
			double dkin= - 2 * kin_i[old_candidate];
			
			dkin+= 2 * (kin_i[new_candidate] - vertices[new_candidate]->links->posweightof(old_candidate).second);
			
			double dkt=vertices[new_candidate]->strength -  vertices[old_candidate]->strength;
			
			double dq=modularity(kin_g+ dkin , ktot_g + dkt) - qtot;
			
			double bol=exp(dq/temp);
	
			if (ran4()<bol) {
				
				kin_g+=dkin;
				ktot_g+=dkt;
				
				group[old_position]=new_candidate;
				
				group_s.erase(old_candidate);
				group_s.insert(new_candidate);
				
				for (int k=0; k<vertices[old_candidate]->links->size(); k++)
					kin_i[vertices[old_candidate]->links->l[k]]-=vertices[old_candidate]->links->w[k];
				
				for (int k=0; k<vertices[new_candidate]->links->size(); k++)
					kin_i[vertices[new_candidate]->links->l[k]]+=vertices[new_candidate]->links->w[k];
					
				
				qtot+=dq;
				
				
				
				
			}
			
			
			
			
						
			
			
		}		
		
		/*
		if (counter%100000==0) {
			cout<<temp<<"\t"<<qtot<<endl;
			cout<<"check .... "<<fabs(qtot - modularity(group))<<endl;
		}
		//*/
		
		if (fabs(qtot-qprevious)< 1e-5)
			stopper++;
		else
			stopper=0;
			
		
		if (stopper==max_loops)
			break;


	}
	
	sort(group.begin(), group.end());
	print_id(group, cout);
	
	
	

	return 0;

}





int clanetwork::good_group_h(const int size, deque<int> & group, double initial_temp) {
	
	
	// the group is built in such a way that it must contain 0, eh eh...
	
	group.clear();
	
	if (size==dim) {
	
		cerr<<"please, it is all the network!"<<endl;
		return -1;
		
	}
	
	good_group_greedy(size, group);
	

	double kin_g=kin(group);
	double ktot_g=ktot(group);
	
	double kin_i[dim];
	for (int i=0; i<dim; i++)
		kin_i[i]=vertices[i]->kplus(group);
	
	set<int> group_s;
	for (int i=0; i<group.size(); i++)
		group_s.insert(group[i]);
		
	
		
	double tstep=0.9999;
	int max_loops=3000;
	
	
	double temp=initial_temp;
	double kprevious=0;
	int stopper=0;
	
	int counter=0;

	while (true) {
		
		temp=temp*tstep;
		kprevious=kin_g;
		
		for (int h=0; h<dim; h++) {
			
			counter++;
			int new_candidate=irand(dim-1);
			while(group_s.find(new_candidate)!=group_s.end()) 
				new_candidate=irand(dim-1);
			
			
			int old_position=irand(group.size()-1);
			int old_candidate=group[old_position];
				
			double dk= - 2 * kin_i[old_candidate];
			
			dk+= 2 * (kin_i[new_candidate] - vertices[new_candidate]->links->posweightof(old_candidate).second);
			
			double bol=exp(dk/temp);
	
			if (ran4()<bol) {
				
				kin_g+=dk;
				group[old_position]=new_candidate;
				
				group_s.erase(old_candidate);
				group_s.insert(new_candidate);
				
				for (int k=0; k<vertices[old_candidate]->links->size(); k++)
					kin_i[vertices[old_candidate]->links->l[k]]-=vertices[old_candidate]->links->w[k];
				
				for (int k=0; k<vertices[new_candidate]->links->size(); k++)
					kin_i[vertices[new_candidate]->links->l[k]]+=vertices[new_candidate]->links->w[k];
			
			}
			
			
			
			
						
			
			
		}		
		
		//if (counter%100000==0)
		//	cout<<temp<<"\t"<<kin_g<<endl;

		
		if (fabs(kin_g-kprevious)< 1e-5)
			stopper++;
		else
			stopper=0;
			
		
		if (stopper==max_loops)
			break;


	}
	
	//cout<<"steps taken= "<<counter<<endl;
	sort(group.begin(), group.end());
	//cout<<"good group"<<endl;
	//prints(group);
	


	return 0;

}




int clanetwork::good_group_greedy(const int size, deque<int> & group) {
	
		
	group.clear();
	
	multimap<int, int> kin_node;
	for (int i=0; i<dim; i++)
		kin_node.insert(make_pair(0, i));
	
	map<int, multimap<int, int>::iterator> node_iterator;
	for (multimap<int, int>::iterator it=kin_node.begin(); it!=kin_node.end(); it++)
		node_iterator.insert(make_pair(it->second, it));
	
	
	int kin_g=0;
	
	int best_node=irand(dim-1);
	while(true) {
		
		
		group.push_back(best_node);
		if (group.size()==size)
			break;
		
		kin_node.erase(node_iterator.find(best_node)->second);
		node_iterator.erase(best_node);
		
		for (int i=0; i<vertices[best_node]->links->size(); i++) {
			
			int neigh=vertices[best_node]->links->l[i];
			
			
			map<int, multimap<int, int>::iterator>::iterator itf = node_iterator.find(neigh);
			
			if (itf!=node_iterator.end()) {
				
				int kin_old=itf->second->first;
				kin_node.erase(itf->second);
				itf->second=kin_node.insert(make_pair(kin_old+cast_int(vertices[best_node]->links->w[i]), neigh));
			
			}
		}
		
		

		deque<int> bests;
		
		multimap<int, int>::iterator itlast=kin_node.end();
		itlast--;
		
		int kin_best = itlast->first;
		
		for (multimap<int, int>::iterator it=itlast; it!=kin_node.begin(); it--) {
		
			if (it->first==kin_best)
				bests.push_back(it->second);
			else
				break;
		}
		
		
		best_node=bests[irand(bests.size()-1)];
		kin_g+= 2 * kin_best;
		//cout<<kin_best<<endl;
		//prints(bests);
		
	}
	
	
	
	
	//cout<<kin_g<<" "<<kin(group)<<endl;
	sort(group.begin(), group.end());
	

	return 0;



}


int clanetwork::just_the_worst(const set<int> & g) {
		
	
	deque<pair<double, int> >  rs;
	compute_rs(g, rs);
	
	
	return rs[rs.size()-1].second;



}



void clanetwork::compute_rs(const set<int> & group, deque<pair<double, int> > & rs) {
	
	rs.clear();
	
	
	double kin_g=kin(group);
	double ktot_g=ktot(group);
	int kout_g= cast_int(ktot_g - kin_g);
	int tm = cast_int(2 * tstrength - ktot_g);
	
	

	for (set<int>::iterator its=group.begin(); its!=group.end(); its++) {
		
		
		int kin_i=cast_int(vertices[*its]->kplus(group));
		int strength_i=cast_int(vertices[*its]->strength);
		
		
		double b1=compute_r_hyper(kin_i, kout_g + 2 * kin_i - strength_i, tm + strength_i, strength_i);
		double b2=compute_r_hyper(kin_i+1, kout_g + 2 * kin_i - strength_i, tm + strength_i, strength_i);
		
		
		double r= b2 + ran4() * (b1 - b2);
		
		rs.push_back(make_pair(r, *its));
		
	}

	sort(rs.begin(), rs.end());
	
	
	
	
	
	

}



double clanetwork::random_cscore(const set<int> & group_s) {
	
	
	//cout<<"size "<<group_ggg.size()<<endl;
	
	if(group_s.size()<3)
		return 1;
	
	
	deque<pair<double, int> >  rs;
	compute_rs(group_s, rs);
	
	//for(int i=0; i<rs.size(); i++)
		//cout<<rs[i].first<<" "<<vertices[rs[i].second]->id_num<<endl;
	
	
	
	double z= rs[rs.size()-1].first;
	double rtilda= rs[rs.size()-2].first;
	//cout<<"z: "<<z<<" rtilda: "<<rtilda<<" rtilda - z: "<<rtilda - z<<endl;

	double fraction=(1-z)/(1-rtilda);
	
	double random_cscore= 1 - pow(fraction, dim - group_s.size() + 1);
	
	
	return random_cscore;


}


double clanetwork::r_score(const set<int> & group_s, double & err) {
	
	
		
	if(group_s.size()<3)
		return 1;
	
	deque<double> averaged;
	for(int i=0; i<100; i++)
		averaged.push_back(random_cscore(group_s));
	
	
	//cout<<"AVERAGED "<<average_func(averaged)<<" +/- "<<sqrt(variance_func(averaged)/(averaged.size()-1))<<endl;
	
	double rscore = average_func(averaged);
	err=sqrt(variance_func(averaged)/(averaged.size()-1));

	return max(0., rscore);


}


double clanetwork::r_score(const deque<int> & group_ggg, double & err) {
	
	
		
	if(group_ggg.size()<3)
		return 1;
	
	
	set<int> group_s;
	
	for (int i=0; i<group_ggg.size(); i++)
		group_s.insert(group_ggg[i]);


	return r_score(group_s, err);


}


double clanetwork::r_score(const deque<int> & group_ggg) {
	
	double err;
	return r_score(group_ggg, err);

}

double clanetwork::r_score(const set<int> & group_ggg) {
	
	double err;
	return r_score(group_ggg, err);

}

//*

void clanetwork::compute_rs_external(const set<int> & group, deque<int> & border_nodes, deque<pair<double, int> > & rs) {
	
	rs.clear();
	
	
	double kin_g=kin(group);
	double ktot_g=ktot(group);
	int kout_g= cast_int(ktot_g - kin_g);
	int tm = cast_int(2 * tstrength - ktot_g);
	
	

	for (deque<int>::iterator its=border_nodes.begin(); its!=border_nodes.end(); its++) {
		
		
		int kin_i=cast_int(vertices[*its]->kplus(group));
		int strength_i=cast_int(vertices[*its]->strength);
		
		
		double b1=compute_r_hyper(kin_i, kout_g, tm, strength_i);
		double b2=compute_r_hyper(kin_i+1, kout_g, tm, strength_i);
		
		double r= b2 + ran4() * (b1 - b2);
		
		rs.push_back(make_pair(r, *its));
		
	}

	sort(rs.begin(), rs.end());
	
	
	

}



/*
double compute_cumulative_sum(int n, int N, double z, double rtilda) {

	// this routine computes the probability cumulative of sum of the bottom n order statistics (N in the sample)
	
	//deque<double> all_realizations;
	
	
	int realizations=10000;
	int good=0;
	
	for(int reals=0; reals<realizations; reals++) {
	
		double sum=0;
		double previous=rtilda;
		for(int i=0; i<n; i++) {
			previous=1. - (1. - previous) * pow(1-ran4(), 1./double(N-i));
			sum+=previous;
			if(sum>z)
				break;
		}
		
		//all_realizations.push_back(sum);
		
		if(sum<z)
			good++;
			
	}
	
	
	//cout<<"Average Sum "<<average_func(all_realizations)<<" +/- "<<sqrt(variance_func(all_realizations))<<endl;
	//cout<<"z & good "<<z<<" "<<good<<endl;
	return double(good)/double(realizations);




}

*/


double clanetwork::compute_score_border(set<int> & group_s, deque<int> & border_nodes) {
	
	
	// this function takes group_s, computes the worst node, insert it into border_nodes.
	// compute all the rs for the border. compare them with rs and compute the score.
	
	
	if(group_s.size()<2)
		return -1;
	
	
	
	deque<pair<double, int> >  rs;
	compute_rs(group_s, rs);
	
	//for(int i=0; i<rs.size(); i++)
		//cout<<rs[i].first<<" "<<vertices[rs[i].second]->id_num<<endl;
	
	int wn= rs[rs.size()-1].second;
	
	
	double z= rs[rs.size()-1].first;
	double rtilda= rs[rs.size()-2].first;
	
	group_s.erase(wn);		// I erased the worst node here!!!!
	border_nodes.push_back(wn);
	
	deque<pair<double, int> >  rsborder;
	compute_rs_external(group_s, border_nodes, rsborder);
	
	
	//cout<<"zeta: "<<z<<" rtilda "<<rtilda<<endl;
	//cout<<"rsborder ********"<<endl;
	
	
	
	if(rsborder.size()>0 && rsborder[0].first<rtilda) {
		
		// here I swapped 
		
		
		double rswap=rsborder[0].first;
		rsborder[0].first=rtilda;
		rtilda=rswap;
		
		//cout<<"swapped new rtilda "<<rtilda<<endl;

	
	}
	
	//for(int i=0; i<rsborder.size(); i++)
		//cout<<rsborder[i].first<<" "<<vertices[rsborder[i].second]->id_num<<endl;

	
	double zsum= 0;
	for(int i=0; i<rsborder.size(); i++)
		zsum+=rsborder[i].first;
	
	//zsum+=z;
	
	
	
	
	
	
	//cout<<"Sumet "<<zsum/border_nodes.size()<<" "<<dim-group_s.size()<<" "<<border_nodes.size()<<" "<<rsborder.size()<<" "<<rtilda<<endl;
	//cerr<<border_nodes.size()<<" "<<rsborder[rsborder.size()-1].first<<endl;
	
	
	
	
	return cumformula_gaussian(zsum, dim - group_s.size(), border_nodes.size(), rtilda);


}


double clanetwork::comp_border(set<int> group_s, double & minbsize) {
	
	
	//ofstream bordout("border_scores.dat");
	
		
	deque<int> border;
	double minbss=2;
	minbsize=group_s.size();
	
	
	while(true) {
	
		double cborder=compute_score_border(group_s, border);
		if(cborder<-0.5)
			break;
				
				
		
		//cout<<border.size()<<" "<<cborder<<endl;
		//bordout<<border.size()<<" "<<cborder<<endl;
		
		if(cborder<minbss) {
			
			minbss=cborder;
			minbsize=border.size();
			
		}
		
		
		//int eeee;
		//cout<<"press to continue"<<endl;
		//cin>>eeee;
	
	
	
	}
	
	//system("cat border_scores.dat");
	//cout<<"MIN "<<minbss<<" "<<minbsize<<endl;
	
	return minbss;




}

double clanetwork::b_score(const set<int> & group_s) {

	double err;
	return b_score(group_s, err);


}

double clanetwork::b_score(const set<int> & group_s, double & err) {
	
	
		
	
	
	if(group_s.size()<3)
		return 1;
	if(double(group_s.size())/dim>0.95)
		cerr<<"WARNING: this module is more than 95% of the whole network-\nour approximations are likely to fail and the b-score might be not reliable"<<endl;
	
	
	
	int reals=20;
	deque<double> minbss;
	deque<double> minbsizes;
	for(int h=0; h<reals; h++) {
		
		double minbsize;
		double minbs=comp_border(group_s, minbsize);
		minbss.push_back(minbs);
		minbsizes.push_back(minbsize);
		
		//cout<<h<<" "<<minbs<<" "<<minbsize<<endl;
	
	
	}

	
	
	double minborderscore=average_func(minbss);
	
	//cout<<"bscore "<<minborderscore<<" +/- "<<sqrt(variance_func(minbss)/minbss.size())<<endl;
	//cout<<"minsize average "<<average_func(minbsizes)<<" +/- "<<sqrt(variance_func(minbsizes)/minbsizes.size())<<endl;
	//cout<<"bscore "<<minborderscore<<" +/- "<<sqrt(variance_func(minbss)/(minbss.size()-1))<<"\tbsize "<<average_func(minbsizes)<<" +/- "<<sqrt(variance_func(minbsizes)/(minbsizes.size()-1))<<endl;

	err=sqrt(variance_func(minbss)/(minbss.size()-1));
	
	
	
	return max(0., minborderscore);


}


double clanetwork::b_score(const deque<int> & group_ggg, double & err) {
	
	
	set<int> group_s;
	
	for (int i=0; i<group_ggg.size(); i++)
		group_s.insert(group_ggg[i]);

	
	
	return b_score(group_s, err);


}




int clanetwork::c__inside(const deque<int> & group_d, double threshold) {
	
	deque<int> erased;
	if(P1.verbose)
		cout<<"looking for a "<<threshold<<" core (using c-score)"<<endl;
	
	set<int> group_s;
	
	for (int i=0; i<group_d.size(); i++)
		group_s.insert(group_d[i]);

	
	
	double cs=c_score(group_s);
	while(cs>threshold && group_s.size()>2) {
		
		int wn=just_the_worst(group_s);
		cs=c_score(group_s);
		if(P1.verbose)
			cout<<"c_inside-> worst node = "<<id_of(wn)<<" cscore "<<cs<<" kin_i/k_i: "<<vertices[wn]->kplus(group_s)<<" / "<<vertices[wn]->strength<<endl;
		if(cs>threshold) {
			
			group_s.erase(wn);
			erased.push_back(wn);
		
		}
		
		
	}
	
	if(erased.size()>0 && group_s.size()>2) {
		
		
		if(P1.verbose)
			cout<<"Erased nodes "<<erased.size()<<endl;
		ofstream outt("c_erased.dat", ios::app);
		print_id(erased, outt);
		if(P1.verbose)
			print_id(erased, cout);

	
	}
	
	if(erased.size()==0) {
		ofstream outt("c_erased.dat", ios::app);
		outt<<endl;
	}
	
	if(group_s.size()>2) {
		
		cout<<"C-"<<P1.threshold*100<<"% core of "<<group_s.size()<<" nodes"<<endl;
		
		ofstream outt("c_core.dat", ios::app);
		print_id(group_s, outt);

		if(P1.verbose)
			print_id(group_s, cout);
		if(P1.verbose)
			cout<<"*******************************************"<<endl;
	
	
	}
	else {
		cout<<"No core found"<<endl;
		ofstream outt("c_erased.dat", ios::app);
		print_id(group_d, outt);
		ofstream outt2("c_core.dat", ios::app);
		outt2<<endl;
		
		return 0;
	}

		
	
	return group_s.size();


}



int clanetwork::r__inside(const deque<int> & group_d, double threshold) {
	
	deque<int> erased;
	if(P1.verbose)
		cout<<"looking for a "<<threshold<<" core (using fast c-score)"<<endl;
	
	set<int> group_s;
	
	for (int i=0; i<group_d.size(); i++)
		group_s.insert(group_d[i]);

	
	
	double rs=r_score(group_s);
	
	while(rs>threshold && group_s.size()>2) {
		
		int wn=just_the_worst(group_s);
		rs=r_score(group_s);
		if(P1.verbose)
			cout<<"c_inside-> worst node = "<<id_of(wn)<<" cscore "<<rs<<" kin_i/k_i: "<<vertices[wn]->kplus(group_s)<<" / "<<vertices[wn]->strength<<endl;		
		if(rs>threshold) {
			
			group_s.erase(wn);
			erased.push_back(wn);
		
		}
		
		
	}
	
	if(erased.size()>0 && group_s.size()>2) {
		
		if(P1.verbose)
			cout<<"Erased nodes "<<erased.size()<<endl;
		if(P1.verbose)
			print_id(erased, cout);
		
		ofstream outt("c_erased.dat", ios::app);
		print_id(erased, outt);


	
	}
	
	if(erased.size()==0) {
		ofstream outt("c_erased.dat", ios::app);
		outt<<endl;
	}
	
	if(group_s.size()>2) {
		
		cout<<"C-"<<P1.threshold*100<<"% core of "<<group_s.size()<<" nodes"<<endl;
		
		ofstream outt("c_core.dat", ios::app);
		print_id(group_s, outt);

		if(P1.verbose)
			print_id(group_s, cout);
		if(P1.verbose)
			cout<<"*******************************************"<<endl;
	
	
	}
	else {
		cout<<"No core found"<<endl;
		ofstream outt("c_erased.dat", ios::app);
		print_id(group_d, outt);
		ofstream outt2("c_core.dat", ios::app);
		outt2<<endl;
		
		return 0;
	}

		
	
	return group_s.size();


}




int clanetwork::b__inside(const deque<int> & group_d, double threshold) {
	
	deque<int> erased;
	if(P1.verbose)
		cout<<"looking for a "<<threshold<<" core (using b-score)"<<endl;
	
	set<int> group_s;
	
	for (int i=0; i<group_d.size(); i++)
		group_s.insert(group_d[i]);

	
	
	double bs=b_score(group_s);
	
	while(bs>threshold && group_s.size()>2) {
		
		int wn=just_the_worst(group_s);
		bs=b_score(group_s);
		if(P1.verbose)
			cout<<"b_inside-> worst node = "<<id_of(wn)<<" bscore "<<bs<<" kin_i/k_i: "<<vertices[wn]->kplus(group_s)<<" / "<<vertices[wn]->strength<<endl;		
		if(bs>threshold) {
			
			group_s.erase(wn);
			erased.push_back(wn);
		
		}
		
		
	}
	
	if(erased.size()>0 && group_s.size()>2) {
		
		if(P1.verbose)
			cout<<"Erased nodes "<<erased.size()<<endl;
		ofstream outt("b_erased.dat", ios::app);
		print_id(erased, outt);
		if(P1.verbose)
			print_id(erased, cout);

	
	}
	
	if(erased.size()==0) {
	
		ofstream outt("b_erased.dat", ios::app);
		outt<<endl;
	}

	if(group_s.size()>2) {
		
		cout<<"B-"<<P1.threshold*100<<"% core of "<<group_s.size()<<" nodes"<<endl;
		
		ofstream outt("b_core.dat", ios::app);
		print_id(group_s, outt);

		if(P1.verbose)
			print_id(group_s, cout);
		if(P1.verbose)
			cout<<"*******************************************"<<endl;
	
	
	}
	else {
	
		cout<<"No core found"<<endl;
		ofstream outt("b_erased.dat", ios::app);
		print_id(group_d, outt);
		ofstream outt2("b_core.dat", ios::app);
		outt2<<endl;
		
		return 0;
	}

		
			
	
	return group_s.size();
	
	
	


}
