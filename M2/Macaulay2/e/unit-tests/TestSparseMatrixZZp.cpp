#include "sparse-matrix/SparseMatrixZZp.hpp"
#include "sparse-matrix/PivotHelper.hpp"
#include "mat-linalg.hpp"
#include "timing.hpp"

#include <iostream>
#include <fstream>
#include <gtest/gtest.h>


TEST(SparseMatrixZZp, fromTriple)
{
  M2:: ARingZZpFlint F(101);
  SparseMatrixZZp A(F, 5, 7, {{0,1,1}, {0,4,2}, {0,6,3},
                           {1,3,4},
                           {2,0,5}, {2,2,6}, {2,6,7},
                           {4,0,8},{4,5,9}});
  A.dump(std::cout);
  
  for (auto r=0; r < A.numRows(); ++r)
    {
      for (auto i = A.cbegin(r); i != A.cend(r); ++i)
        {
          auto t = *i;
          std::cout << "[" << r << "," << t.first << "," << t.second << "] ";
        }
      std::cout << std::endl;
    }

  A.denseDisplay(std::cout);
}

TEST(SparseMatrixZZp, trivialPivotTest)
{
  M2:: ARingZZpFlint F(101);
  SparseMatrixZZp A(F, 5, 7, {{ 0,1,1}, {0,4,2}, {0,6,3},
                           {1,3,4},
                           {2,0,5}, {2,2,6}, {2,6,7},
                           {4,0,8},{4,5,9}});
  PivotHelper pivotHelper(A);
  pivotHelper.findTrivialRowPivots(A);
  std::cout << pivotHelper << std::endl;
  pivotHelper.findTrivialColumnPivots(A);
  std::cout << pivotHelper << std::endl;

  SparseMatrixZZp A2(F, 6, 9, {{0,0,1}, {0,2,2}, {0,5,3}, {0,6,4}, {0,8,5},
  			      {1,1,1}, {1,2,2}, {1,3,3}, {1,5,4}, {1,7,5},
  			      {2,2,1}, {2,3,2}, {2,7,3}, {2,8,4},
  			      {3,1,1}, {3,2,2}, {3,4,3},
  			      {4,1,1}, {4,3,2}, {4,6,3}, {4,8,4},
  			      {5,0,1}, {5,2,2}, {5,4,3}, {5,5,4}, {5,7,5}, {5,8,6}});
  A2.denseDisplay(std::cout);
  
  PivotHelper pivotHelper2(A2);
  pivotHelper2.findTrivialRowPivots(A2);
  std::cout << pivotHelper2 << std::endl;
  pivotHelper2.findTrivialColumnPivots(A2);
  std::cout << pivotHelper2 << std::endl;
  pivotHelper2.findRemainingPivotsGreedy(A2);
  std::cout << pivotHelper2 << std::endl;

  std::vector<IndexType> rowPerm;
  std::vector<IndexType> columnPerm;
  
  pivotHelper2.findUpperTrapezoidalPermutations(A2,rowPerm,columnPerm);

  SparseMatrixZZp B = A2.applyPermutations(rowPerm,columnPerm);
  B.denseDisplay(std::cout);

}

#if 0

TEST(SparseMatrixZZp, fromTriplesFile)
{
  M2:: ARingZZpFlint F(101);
  std::ifstream infile;
  infile.open("sparsemat-1.txt");
  if (not infile)
    {
      std::cout << "file not opened properly" << std::endl;
      exit(1);
    }
  
  SparseMatrixZZp A(F, infile);

  infile.close();;

  SparseMatrixZZp B(F, 5, 7, {{0,1,1}, {0,4,2}, {0,6,3},
                           {1,3,4},
                           {2,0,5}, {2,2,6}, {2,6,7},
                           {4,0,8}});

  A.denseDisplay(std::cout);
  B.denseDisplay(std::cout);

  // eventually: check that A == B.
}

//#endif

TEST(SparseMatrixZZp, fromTriplesFile2)
{
  M2::ARingZZpFlint F(101);
  M2::ARingZZpFFPACK F_ffpack(101);
  std::ifstream infile;
  infile.open("../exampleMat.sms");
  if (not infile)
    {
      std::cout << "file not opened properly" << std::endl;
      exit(1);
    }
  
  SparseMatrixZZp A(F, infile);
  //A.denseDisplay(std::cout);

  auto tt = now();
  DMat<M2::ARingZZpFFPACK> A_DMat(F_ffpack,A.numRows(),A.numColumns());
  std::cout << "Time spent creating DMat: " << seconds(now() - tt) << std::endl;
  toDenseMatrix(A,F,A_DMat);
  auto a = DMatLinAlg<M2::ARingZZpFFPACK>(A_DMat);
  tt = now();
  long myRank = a.rank();
  std::cout << "Time spent computing rank: " << seconds(now() - tt) << std::endl;
  std::cout << "Rank of matrix   : " << myRank << std::endl;


  std::cout << "Number of rows   : " << A.numRows() << std::endl;
  std::cout << "Number of columns: " << A.numColumns() << std::endl;

  PivotHelper pivotHelper(A);

  auto t1 = now();
  pivotHelper.findPivots(A);
  std::cout << "Time spent finding pivots: " << seconds(now() - t1) << std::endl;
  std::cout << "Number of pivots found: " << pivotHelper.numPivots() << std::endl;

  std::cout << pivotHelper << std::endl;

  infile.close();;

  //  A.denseDisplay(std::cout);
}

#endif
//#if 0

TEST(SparseMatrixZZp, fromTriplesFile3)
{
  M2::ARingZZpFlint F(101);
  M2::ARingZZpFFPACK F_ffpack(101);
  std::ifstream infile;
  //infile.open("/Users/moorewf/Downloads/10164x1740.sms");
  //infile.open("/Users/moorewf/Downloads/47104x30144bis.sms");
  //infile.open("/Users/moorewf/Downloads/GL7d16.sms");
  infile.open("/Users/frank/Downloads/GL7d15.sms");
  //infile.open("/Users/moorewf/Downloads/GL7d13.sms");
  //infile.open("/Users/moorewf/Downloads/GL7d12.sms");
  //infile.open("../exampleMat.sms");
  if (not infile)
  {
    std::cout << "file not opened properly" << std::endl;
    exit(1);
  }
  
  SparseMatrixZZp B(F, infile);
  infile.close();
  
  std::vector<IndexType> rowPerm(B.numRows());
  std::vector<IndexType> columnPermInverse(B.numColumns());
  
  std::iota(rowPerm.begin(),rowPerm.end(),0);
  for (int i = 0; i < B.numColumns(); ++i)
     columnPermInverse[i] = B.numColumns()-1-i;
  
  SparseMatrixZZp A = B.applyPermutations(rowPerm,columnPermInverse);

  rowPerm.resize(0);
  columnPermInverse.resize(0);

  //SparseMatrixZZp A = B.transpose();

  // auto tt = now();
  // DMat<M2::ARingZZpFFPACK> A_DMat(F_ffpack,A.numRows(),A.numColumns());
  // std::cout << "Time spent creating DMat: " << seconds(now() - tt) << std::endl;
  // toDenseMatrix(A,F,A_DMat);
  // auto a = DMatLinAlg<M2::ARingZZpFFPACK>(A_DMat);
  // tt = now();
  // long myRank = a6.rank();
  // std::cout << "Time spent computing rank: " << seconds(now() - tt) << std::endl;
  // std::cout << "Rank of matrix   : " << myRank << std::endl;

  std::cout << "Number of rows   : " << A.numRows() << std::endl;
  std::cout << "Number of columns: " << A.numColumns() << std::endl;

  PivotHelper pivotHelper(A);

  auto t1 = now();
  pivotHelper.findPivots(A);
  std::cout << "Time spent finding pivots: " << seconds(now() - t1) << std::endl;
  std::cout << "Number of pivots found: " << pivotHelper.numPivots() << std::endl;

  std::cout << "Finding UT Permutations." << std::endl << std::flush;
  t1 = now();
  pivotHelper.findUpperTrapezoidalPermutations(A,rowPerm,columnPermInverse);
  std::cout << "Time spent finding permutations: " << seconds(now()-t1) << std::endl << std::flush;

  std::cout << "Checking UT Permutations... ";
  if (!A.checkUpperTrapeziodalPermutations(rowPerm,columnPermInverse,pivotHelper.numPivots()))
     std::cout << "Error in finding pivots and/or permutations." << std::endl;
  else
     std::cout << "Passed!" << std::endl;
  //  A.denseDisplay(std::cout);
}

//      if q: 4 2 1 6 5 0 3

//            0 1 2 3 4 5 6

//then qinv : 5 2 1 6 0 4 3

//#endif

TEST(SparseMatrixZZp, fromUnsortedTriple)
{
  M2:: ARingZZpFlint F(101);
  SparseMatrixZZp A(F, 5, 7, {{0,6,3},
                           {1,3,4}, {2,6,7},
                           {2,0,5}, {2,2,6}, {0,4,2}, 
                           {0,1,1}, {4,0,8}});
  A.dump(std::cout);
  
  for (auto r=0; r < A.numRows(); ++r)
    {
      for (auto i = A.cbegin(r); i != A.cend(r); ++i)
        {
          auto t = *i;
          std::cout << "[" << r << "," << t.first << "," << t.second << "] ";
        }
      std::cout << std::endl;
    }

  A.denseDisplay(std::cout);
}

TEST(SparseMatrixZZp, zeroMatrix)
{
  M2:: ARingZZpFlint F(101);
  SparseMatrixZZp A(F, 0, 0, {{}});

  for (auto r=0; r < A.numRows(); ++r)
    {
      for (auto i = A.cbegin(r); i != A.cend(r); ++i)
        {
          auto t = *i;
          std::cout << "[" << r << "," << t.first << "," << t.second << "] ";
        }
      std::cout << std::endl;
    }
}

TEST(SparseMatrixZZp, zeroEntries)
{
  M2:: ARingZZpFlint F(101);
  SparseMatrixZZp A(F, 3, 5, {});

  A.dump(std::cout);
  
  for (auto r=0; r < A.numRows(); ++r)
    {
      for (auto i = A.cbegin(r); i != A.cend(r); ++i)
        {
          auto t = *i;
          std::cout << "[" << r << "," << t.first << "," << t.second << "] ";
        }
      std::cout << std::endl;
    }
}

TEST(SparseMatrixZZp, transpose)
{
  M2:: ARingZZpFlint F(101);
  SparseMatrixZZp A(F, 5, 7, {{0,1,1}, {0,4,2}, {0,6,3},
                           {1,3,4},
                           {2,0,5}, {2,2,6}, {2,6,7},
                           {4,0,8}});

  A.denseDisplay(std::cout);

  std::cout << "Now for the transpose:" << std::endl;
  SparseMatrixZZp B = A.transpose();
  B.denseDisplay(std::cout);
}

TEST(SparseMatrixZZp, randomSparseMatrix)
{
  M2:: ARingZZpFlint F(101);
  auto A = SparseMatrixZZp::randomSparseMatrix(F, 50, 60, .1);
  // A.denseDisplay(std::cout);
  auto B = A.transpose();
  // std::cout << "and its transpose is: " << std::endl;
  // B.denseDisplay(std::cout);
}

// TEST(SparseMatrixZZp, bigRandomSparseMatrix)
// {
//   M2:: ARingZZpFlint F(101);
//   auto A = SparseMatrixZZp::randomSparseMatrix(F, 50000, 60000, .01);
//   std::cout << "and now lets work on the transpose" << std::endl;
//   auto t1 = now();
//   auto B = A.transpose();
//   std::cout << "transose time: " << seconds(now() - t1) << std::endl;
// }
