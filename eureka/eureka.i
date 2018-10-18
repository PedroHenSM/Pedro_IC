/* eureka.i */
%module eureka
%include "carrays.i" /* %array_functions(type,name) */
%array_functions (double, doubleArray);
/* Creates an void array | array_functions(type,name) */
/* %array_class(void,voidArray) Creates an void array for classes | array_class(type,name) */
/* %typemap(in) void* = double*; */
%{
/* Inlcude headers files or function declarations */
#include "EurekaOptimaException.h"
#include "Problem.h"
#include "TrussBarStructureStaticProblem.h"
#include "TrussBarStructureStaticSimulator.h"
#include "F101Truss10Bar.h"
#include "F103Truss25Bar.h"
#include "F105Truss60Bar.h"
#include "F107Truss72Bar.h"
#include "F109Truss942Bar.h"

double **new_doubleddArray(int rows){
	double **arr = new double *[rows];
	return arr;
}

double **castToDouble(void *b){
	return (double**)b;
}


double **new_doubleddArray(int rows, int cols) {
    int i;
    double **arr = new double *[rows];
    for (i=0; i<rows; i++)
		arr[i] = new double[cols];
    return arr;
}

void delete_doubleddArray (double **arr, int rows, int cols){
	int i;
	for (i=0; i<rows; i++)
		delete[] arr[i];
	delete[] arr;
}

void doubleddArray_setitem(double **array, int row, int col, double value) {
    array[row][col] = value;
}

double doubleddArray_getitem(double **array, int row, int col) {
    return array[row][col];
}
%}

// Helper function to create a 2d array

double **new_doubleddArray(int rows);

double **castToDouble(void *b);

void delete_doubleddArray (double **arr, int rows, int cols);

void doubleddArray_setitem(double **array, int row, int col, double value);

double doubleddArray_getitem(double **array, int row, int col);

/* Inserts class */


namespace problem {

    class Problem {
    public:
        Problem(int dimension, void* bounds, int maxNumberObjectiveFunctionEvaluations, int numberObjectives, int numberConstraints);
        Problem(const Problem& orig);
        virtual ~Problem();
        void evaluate(void* vector, void* values);
		// void evaluate(void** vectors, unsigned int vectorsLength, void* values, unsigned int valuesLength);
        int getNumberObjectiveFunctionEvaluations() const;
        int getMaxNumberObjectiveFunctionEvaluations() const;
        void* getBounds() const;
        int getNumberObjectives() const;
        int getNumberConstraints() const;
        int getDimension() const;
		virtual string toString() const;
    protected:
        int dimension;
        virtual void evaluation(void* vector, void* values) = 0;
		//TODO - Should I include this here or in the ParallelEurekaOptima project ?
		//virtual void evaluation(void** vectors, unsigned int vectorsLength, void** values, unsigned int valuesLength) = 0;
        int numberObjectiveFunctionEvaluations;
        int maxNumberObjectiveFunctionEvaluations;
        void* bounds;
        int numberObjectives;
        int numberConstraints;
    };

    // typedef shared_ptr<Problem> ProblemPtr; // TODO: MODIFIED. 'shared_ptr' is ambiguous
    typedef std::shared_ptr<Problem> ProblemPtr;
}

namespace problem {

	class TrussBarStructureStaticProblem: public Problem {
	public:
		TrussBarStructureStaticProblem(int dimension, void* bounds, int maxNumberObjectiveFunctionEvaluations, int numberObjectives, int numberConstraints, int numberOfBars, double gamma, double stressConstraint, double displacementConstraint, /*int numberOfDisplacementConstraints,*/ string inputFileName, double lowerBound, double upperBound);
		TrussBarStructureStaticProblem(const TrussBarStructureStaticProblem& orig);
		virtual ~TrussBarStructureStaticProblem();
		virtual void evaluation(void* vector, void* values);
		int const getNumberOfBars();
		int const getNLCase();
		double const getDisplacementConstraint();
	protected:
		/*
		 * Return an integer array containing information about how to group the bars of the structure.
		 * NULL is indicated to the cases in which there is not grouping, that is, when the number of bars is equals to the number of groups.
		 */
		virtual int* const getGrouping();
		/*
		 * Create a array with the areas of the bars using the grouping.
		 * It is important highlight that the method 'createGrouping' should be used here.
		 */
		virtual void fillAreasAux(double* x);
		/*
		 * Return the array with the cross-section areas.
		 */
		double* const getAreasAux();
		double* const getStressDisplacementAux();
		/*
		 * Calculate the number of displacements.
		 */
		int getNumberOfDisplacements();
		/*
		 * Returns a pointer to the simulator of truss structures.
		 */
		TrussBarStructureStaticSimulatorPtr const getSimulator();
	private:
		/*
		 * Simulator used to calculate the weight, stresses and displacements of the structure.
		 */
		TrussBarStructureStaticSimulatorPtr simulator;
		/*
		 * Number of bars of the structure.
		 */
		int numberOfBars;
		/*
		 * Specific weight of the bar's material.
		 */
		double gamma;
		/*
		 * Maximum stress allowed in the bars.
		 */
		double stressConstraint;
		/*
		 * Maximum displacement allowed for the nodes.
		 */
		double displacementConstraint;
		/*
		 * Number of nodes where displacements can happen.
		 */
		int numberOfDisplacements;
		/*
		 * Variable to archive the stresses and displacements.
		 */
		double* stressDisplacementAux;
		/*
		 * Variable to archive the cross-sectional area of the bars.
		 */
		double* areasAux;
		//variables with data from input file
		int numnp_cpp;
		int numeg_cpp;
		int nlcase_cpp;
		int modex_cpp;
		int* node_n_cpp;
		int** node_id_cpp;
		double** node_position_cpp;
		int* ll_cpp;
		int* loads_n_cpp;
		int** loads_node_cpp;
		int** loads_direction_cpp;
		double** loads_cpp;
		int elements_npar_1;
		int elements_npar_2;
		int elements_npar_3;
		int* element_id_cpp;
		double* element_cpp;
		double* element_extra_cpp; //this information is not used by stap
		int* m_cpp;
		int* ii_cpp;
		int* jj_cpp;
		int* mtyp_cpp;
		int* kg_cpp;

		/*
		 * Load information about the structure from a file.
		 */
		void readFile(string fileName);

	};

}

namespace problem {

    class F101Truss10Bar: public TrussBarStructureStaticProblem {
    public:
        F101Truss10Bar();
        F101Truss10Bar(const F101Truss10Bar& orig);
        virtual ~F101Truss10Bar();
    };

}

namespace problem {

    class F103Truss25Bar: public TrussBarStructureStaticProblem {
    public:
        F103Truss25Bar();
        F103Truss25Bar(const F103Truss25Bar& orig);
        virtual ~F103Truss25Bar();
		virtual void evaluation(void* vector, void* values);
	protected:
		virtual int* const getGrouping();
    private:
        int* grouping;
    };

}

namespace problem {

    class F105Truss60Bar: public TrussBarStructureStaticProblem {
    public:
        F105Truss60Bar();
        F105Truss60Bar(const F105Truss60Bar& orig);
        virtual ~F105Truss60Bar();
		virtual void evaluation(void* vector, void* values);
	protected:
        virtual int* const getGrouping();
		virtual void fillAreasAux(double* x);
    private:
        int* grouping;
		double displacementConstraint2;
		double displacementConstraint3;
    };

}

namespace problem {

    class F107Truss72Bar: public TrussBarStructureStaticProblem {
    public:
        F107Truss72Bar();
        F107Truss72Bar(const F107Truss72Bar& orig);
        virtual ~F107Truss72Bar();
		virtual void evaluation(void* vector, void* values);
	protected:
        virtual int* const getGrouping();
    private:
        int* grouping;
    };

}

namespace problem {

    class F109Truss942Bar: public TrussBarStructureStaticProblem {
    public:
        F109Truss942Bar();
        F109Truss942Bar(const F109Truss942Bar& orig);
        virtual ~F109Truss942Bar();
    protected:
        virtual int* const getGrouping();
		virtual void fillAreasAux(double* x);
	private:
		int* grouping;
    };

}
