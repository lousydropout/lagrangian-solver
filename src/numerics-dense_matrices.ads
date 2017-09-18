with Numerics;
use  Numerics;
package Numerics.Dense_Matrices is
   
   Error : exception;
   
   function Outer (X, Y : in Real_Vector) return Real_Matrix;
   function "*" (X, Y : in Real_Vector) return Real_Matrix renames Outer;
   
   function "*" (X : in Real;
		 A : in Real_Matrix) return Real_Matrix;
   function "*" (A : in Real_Matrix;
		 X : in Real) return Real_Matrix is (X * A);
   
   function "+" (A : in Real_Matrix) return Real_Matrix is (A);
   function "-" (A : in Real_Matrix) return Real_Matrix;
   
   function "+" (A : in Real_Matrix;
		 B : in Real_Matrix) return Real_Matrix
     with Pre => A'Length (1) = B'Length (1) and A'Length (2) = B'Length (2);
     
   function "-" (A : in Real_Matrix;
		 B : in Real_Matrix) return Real_Matrix
     with Pre => A'Length (1) = B'Length (1) and A'Length (2) = B'Length (2);
     
   function "*" (A : in Real_Matrix;
		 B : in Real_Matrix) return Real_Matrix
     with Pre => A'Length (2) = B'Length (1);
   
   function Transpose (A : in Real_Matrix) return Real_Matrix;
   
   function Eye (N : in Pos) return Real_Matrix;
   function Eye (N : in Pos) return Int_Matrix;
   
   procedure LU_Decomposition (A : in     Real_Matrix;
			       P :    out Int_Array;
			       L :    out Real_Matrix;
			       U :    out Real_Matrix)
     with Pre => A'Length (1) = A'Length (2);
   
   procedure Print (X : in Real_Vector);
   procedure Print (A : in Real_Matrix);
   procedure Print (A : in Int_Matrix);
   
   function Determinant (P : in Int_Array;
			 L : in Real_Matrix;
			 U : in Real_Matrix) return Real
     with Pre => L'Length (1) = L'Length (2) 
     and U'Length (1) = U'Length (2) 
     and L'Length (1) = U'Length (1)
     and L'Length (1) = P'Length;
   
   function Determinant (A : in Real_Matrix) return Real
     with Pre => A'Length (1) = A'Length (2);
   
   function Solve_Upper (U : in Real_Matrix;
   			 B : in Real_Vector) return Real_Vector
     with Pre => U'Length (1) = U'Length (2) and U'Length (1) = B'Length;
   
   function Solve_Lower (L : in Real_Matrix;
   			 B : in Real_Vector) return Real_Vector
     with Pre => L'Length (1) = L'Length (2) and L'Length (1) = B'Length;
   
   function Solve (P : in Int_Array;
   		   L : in Real_Matrix;
   		   U : in Real_Matrix;
   		   B : in Real_Vector) return Real_Vector
     with Pre => L'Length (1) = L'Length (2) 
     and U'Length (1) = U'Length (2) 
     and L'Length (1) = U'Length (1)
     and L'Length (1) = P'Length
     and L'Length (1) = B'Length;
   
   function Solve (A : in Real_Matrix;
		   B : in Real_Vector) return Real_Vector
     with Pre => A'Length (1) = A'Length (2)
     and A'Length (1) = B'Length;
   
   function Inverse (A : in Real_Matrix) return Real_Matrix
     with Pre => A'Length (1) = A'Length (2);
   
   function Diag (A : in Real_Matrix) return Real_Vector
     with Pre => A'Length (1) = A'Length (2);
   function Diag (X : in Real_Vector) return Real_Matrix;
   
   function "and" (A, B : in Real_Matrix) return Real_Matrix;
   function "or" (A, B : in Real_Matrix) return Real_Matrix;
   
   function Remove_1st_N (X : in Real_Vector;
			  N : in Pos) return Real_Vector
     with Pre => X'Length > N;
   function Remove_1st_N (A : in Real_Matrix;
			  N : in Pos) return Real_Matrix
     with Pre => A'Length (1) > N and A'Length (2) > N;
   
   procedure Copy (From	 : in     Real_Vector;
		   To	 : in out Real_Vector;
		   Start : in     Pos)
     with Pre => To'Length >= From'Length and To'First <= Start;
   
   procedure Copy (From	   : in     Real_Matrix;
		   To	   : in out Real_Matrix;
		   Start_I : in     Pos;
		   Start_J : in     Pos)
     with Pre => To'Length (1) >= From'Length (1)
     and To'Length (2) >= From'Length (2)
     and To'First (1) <= Start_I
     and To'First (2) <= Start_J;
   
private
   
   function Pivoting_Array (A : in Real_Matrix) return Int_Array;
   function Permute_Row (A : in Real_Matrix;
			 P : in Int_Array) return Real_Matrix
     with Pre => A'Length (1) = P'Length;
   
   function Permute_Row (X : in Real_Vector;
			 P : in Int_Array) return Real_Vector
     with Pre => X'Length = P'Length;
   function Number_Of_Swaps (P : in Int_Array) return Pos;
     
   --  procedure Permute_Col (A : in out Real_Matrix;
   --  			  P : in     Int_Array)
   --    with Pre => A'Length (2) = P'Length;
   
end Numerics.Dense_Matrices;
