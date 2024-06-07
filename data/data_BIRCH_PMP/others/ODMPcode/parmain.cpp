#include "mpi.h"
#include <stdio.h>


#include "logerr.h"
#include "Convert.h"
#include "PARAMFILE.h"
#include "ODMProblem.h"
int main(int argc, char* argv[])
{

	int myid, size;
    int  namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
	
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(processor_name,&namelen);
	
	String s;
    sprintf(s, "Process %d on %s\n", myid, processor_name);
	print(s);


	int parmethod[2];

	if(myid == 0)
	{
		// I am master


		// Cheking program parameters
		//if(argc!=2)
		//{
		//	print("\nConverting files...\n");
		//	error("Enter exactly four program paramenters: parameter file name,\n", 
		//		"Log is writen in .log.");
		//}

		PARAMFILE par(argv[1]);

		String fn1, fn2, fn3, fn4, fn5, fn6;

		par.strParam("probname", fn1);
		par.strParam("namefile", fn2);
		par.strParam("costfile", fn3);
		par.strParam("modelfile", fn4);
		par.strParam("modelnamefile", fn5);
		par.strParam("constrfile", fn6);

		int dp = par.intParam("dp");

		//Initialize object convert, where argv[1] is the problem name
		Convert* convert=new Convert(fn1, fn2, fn3, fn4, fn5, fn6);


		// processing
		convert->process();

		int** solutions=new int*[dp+1];
		int i;
		for(i=0; i<dp+1; i++)
		{
			solutions[i]=new int[convert->nmodels];
		}
		double* values=new double[dp+1];

		//Get data to solve
		ODMProblem* prob= new ODMProblem(convert->probname, convert->nmodels, convert->modelmasks,
										 convert->modelcosts, convert->modeldemands, 
										 convert->realmodelnames, convert->noptions,
										 convert->hardoptions, dp, solutions, values,
										 convert->constructioncost);





		//optimizing
		int method = atoi(argv[2]);//par.intParam("method");


		

		switch(method)
		{
		case 0:
			break;
		case 1:
			parmethod[0] = 0;
			parmethod[1] = 0;
			break;
		case 2:
			parmethod[0] = 1;
			parmethod[1] = 0;
			break;
		case 3:
			parmethod[0] = 1;
			parmethod[1] = 1;
			break;
		case 4:
			parmethod[0] = 1;
			parmethod[1] = 2;
			break;
		default:
			print("Unknown method!");
			break;
		}

		MPI_Bcast(&parmethod, 2, MPI_INT, 0, MPI_COMM_WORLD);

		//String s;
		//sprintf(s, "id=%d lh=%d level=%d\n", myid, parmethod[0], parmethod[1]);
		//print(s);

		prob->optimizeParInit(parmethod[0], parmethod[1]);

		
		// sending/recieving subproblems

		/*for( i=0; i<prob->ngraphcomponents; i++)
		{
			Buffer buf;
			prob->subgraphs[i]->packprob(buf);
			GraphComponent* sg = new GraphComponent(buf);
			delete[] buf.data;

			sg->greedy(parmethod[0], parmethod[1]);

			Buffer buf1;
			sg->packsol(buf1);
			prob->subgraphs[i]->unpacksol(buf1);
			delete[] buf1.data;

		}*/

		
		time_t beg_t = time(NULL);

		int nsent = 0;
		
		int nng = prob->ngraphcomponents;
		int* ns = new int[nng];
		
		int k;
		for(k=0; k < nng; k++)
		{
			ns[k] = prob->subgraphs[k]->nmodels;
			printf("\nns[%3.3d]=%5d\n", k, ns[k]);
		}
		
		for( k=1; k<size; k++)
		{	
			if(nsent == prob->ngraphcomponents)
			{
				MPI_Send(0, 0, MPI_CHAR, k, 0, MPI_COMM_WORLD);
				continue;
			}
			int ind = -1;
			int mind = -1;
			
			int fk;
			for(fk=0; fk < nng; fk++)
			{
				if( mind < ns[fk])
				{
					mind = ns[fk];
					ind = fk;
				}
			}
			ns[ind] = -1;
			//ind = nsent;
			
			printf("\n\nind=%d\n", ind);
			
			Buffer buf;
			prob->subgraphs[ind]->packprob(buf);						
			MPI_Send(buf.data, buf.size, MPI_CHAR, k, buf.size, MPI_COMM_WORLD);
						
			delete[] buf.data;
			nsent++;
		}

		for(k=0; k<prob->ngraphcomponents; k++)
		{
			MPI_Status status;
			MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			//sprintf(s, "source=%d cint=%d size=%d\n", status.MPI_SOURCE, status.MPI_TAG, prob->subgraphs[status.MPI_TAG]->solsize);
			//print(s);
			int cind = status.MPI_TAG;
			Buffer buf;
			buf.size = prob->subgraphs[cind]->solsize;
			//sprintf(s, "size=%d\n", buf.size);
			print(s);
			buf.data = new char[buf.size];

			MPI_Recv(buf.data, buf.size, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
			//print("received\n");
			prob->subgraphs[cind]->unpacksol(buf);
			//print("unpacked\n");
			delete[] buf.data;

			//print("deleted\n");
			int source = status.MPI_SOURCE;
			//sprintf(s, "got %d\n", source);
			print(s);
			if(nsent == prob->ngraphcomponents)
			{
				//sprintf(s, "send to finish %d\n", source);
				//print(s);
				MPI_Send(0, 0, MPI_CHAR, source, 0, MPI_COMM_WORLD);
			}
			else
			{
				//sprintf(s, "send %d\n", source);
				//print(s);
				int ind = -1;
				int mind = -1;
				
				int fk;
				for(fk=0; fk < nng; fk++)
				{
					if( mind < ns[fk])
					{
						mind = ns[fk];
						ind = fk;
					}
				}
				ns[ind] = -1;
				//ind = nsent;
				Buffer buf1;
				prob->subgraphs[ind]->packprob(buf1);						
				MPI_Send(buf1.data, buf1.size, MPI_CHAR, source, buf1.size, MPI_COMM_WORLD);
				delete[] buf1.data;
				nsent++;
			}

		}

		
		prob->optimizeParFinal(parmethod[0], parmethod[1]);
		
		delete[] ns;
		
		time_t end_t = time(NULL);
		printf("*********************************************************\n");
		printf("*********************************************************\n");
		printf("wall time = %d\n", end_t-beg_t);
		printf("*********************************************************\n");
		printf("*********************************************************\n");
		
		std::ofstream to("globrez.log", std::ios::out | std::ios::app);
		
		to << convert->probname << "\t" << "v1" << "\t" << method << "\t" << size-1 << "\t" << end_t-beg_t << "\n";
		
		to.close();

		// saving solution
		//prob->writesol();

		// saving graph
		//prob->writegraph();
		

		delete[] values;
		for( i=0; i<dp+1; i++)
		{
			delete[] solutions[i];
		}
		delete[] solutions;

		//deleting prob
		delete prob;

		//deleting convert
		delete convert;

	}
	else
	{
		// I am slave
		MPI_Bcast(&parmethod, 2, MPI_INT, 0, MPI_COMM_WORLD);
		String s;
		//sprintf(s, "process=%d lh=%d level=%d\n", myid, parmethod[0], parmethod[1]);
		//print(s);

		while(true)
		{
			MPI_Status status;
			MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			if(status.MPI_TAG == 0)
			{
				sprintf(s, "\nprocess=%d finished\n", myid);
				print(s);
				break;
			}
			Buffer buf;
			buf.size = status.MPI_TAG;
			buf.data = new char[buf.size];
			
			MPI_Recv(buf.data, buf.size, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			printf("\n\nProcess %d: ", myid);
			GraphComponent* sg = new GraphComponent(buf);
			delete[] buf.data;

			Buffer buf1;
			//sg->printg();
			sg->greedy(parmethod[0], parmethod[1]);

			sg->packsol(buf1);

			MPI_Send(buf1.data, buf1.size, MPI_CHAR, 0, sg->cind, MPI_COMM_WORLD);
			delete[] buf1.data;

			delete sg;
		}
	}
		
    MPI_Finalize();
	return 0;
}
