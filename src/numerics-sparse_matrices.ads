with Numerics;
use  Numerics;
package Numerics.Sparse_Matrices is
   
   
   -------- Define Enumeration types --------------------------------
   type Permute_By_Type is (Row, Column);
   type Sparse_Matrix_Format   is (CSR, CSC, Triplet);
   package Sparse_Matrix_Format_IO is new Ada.Text_IO.Enumeration_IO (Sparse_Matrix_Format);
   
   ------- Define Matrix --------------------------------------------
   type Sparse_Matrix is private;
   
   
   --- Print procedure ----------------------------------------------
   procedure Print (Mat : in Sparse_Matrix); 
   
   ------- Basic Getter Functions -----------------------------------
   function Norm2 (Item : in Sparse_Matrix) return Real;
   function N_Row (Mat : in Sparse_Matrix)  return Pos;
   function N_Col (Mat : in Sparse_Matrix)  return Pos;
   function Number_Of_Elements (X : in Sparse_Matrix) return Integer;
   
   ------- Functions for Creating Sparse Matrices -------------------
   function Add_Column (X : in Sparse_Matrix;
			V : in Sparse_Vector) return Sparse_Matrix;
   function As_Matrix (X : in Sparse_Vector) return Sparse_Matrix;
   
   function "and" (A, B : in Sparse_Vector) return Sparse_Matrix
     is (Add_Column (As_Matrix (A), B));
   
   function "and" (X : in Sparse_Matrix; V : in Sparse_Vector) 
		  return Sparse_Matrix renames Add_Column;
		   
   procedure Set_Diag (X  : in out Sparse_Matrix;
   		       To : in     Sparse_Vector)
     with Pre => Is_Square_Matrix (X) and N_Col (X) = Length (To);
   function Diag (X : in Sparse_Matrix) return Sparse_Vector
     with Pre => Is_Square_Matrix (X);
   function Diag (X : in Sparse_Vector) return Sparse_Matrix;
   function Sparse (X	: in Real_Matrix;
		    Eps	: in Real	 := 10.0 * Real'Small) 
		   return Sparse_Matrix;
   function Triplet_To_Matrix (I      : in Int_Array;
			       J      : in Int_Array;
			       X      : in Real_Vector;
			       N_Row  : in Pos := 0;
			       N_Col  : in Pos := 0;
			       Format : in Sparse_Matrix_Format := CSC) 
			      return Sparse_Matrix
     with Pre => I'Length = J'Length and I'Length = X'Length;
   function Convert (Mat : in Sparse_Matrix) return Sparse_Matrix;
   procedure Add (Mat  : in out Sparse_Matrix;
		  I, J : in     Nat;
		  X    : in     Real)
     with Pre => I <= N_Row (Mat) and J <= N_Col (Mat);
   procedure Set (Mat  : in out Sparse_Matrix;
		  I, J : in     Nat;
		  X    : in     Real)
     with Pre => I <= N_Row (Mat) and J <= N_Col (Mat);
   
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Testing Functions ----------------------------------------
   function Is_Square_Matrix (A : in Sparse_Matrix) return Boolean;
   function Has_Same_Dimensions (Left, Right : in Sparse_Matrix) return Boolean;
   function Is_Valid (Mat : in Sparse_Matrix) return Boolean;
   
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Matrix operations ----------------------------------------
   function Eye (N : in Nat) return Sparse_Matrix;
   function Zero (N : in Nat) return Sparse_Matrix;
   function Omega (N : in Nat;
		   M : in Pos := 0) return Sparse_Matrix;
   
   function Transpose (Mat : in Sparse_Matrix) return Sparse_Matrix;
   function Plus (Left  : in Sparse_Matrix;
		  Right : in Sparse_Matrix) return Sparse_Matrix
     with Pre => Has_Same_Dimensions (Left, Right);
   function Plus2 (Left  : in Sparse_Matrix;
		   Right : in Sparse_Matrix) return Sparse_Matrix
     with Pre => Has_Same_Dimensions (Left, Right);
   function Minus (Left  : in Sparse_Matrix;
		   Right : in Sparse_Matrix) return Sparse_Matrix
     with Pre => Has_Same_Dimensions (Left, Right);
   function Mult (Left, Right : in Sparse_Matrix) return Sparse_Matrix
     with Pre => N_Col (Left) = N_Row (Right);
   function Kronecker (A, B : in Sparse_Matrix) return Sparse_Matrix;
   function Direct_Sum (A, B : in Sparse_Matrix) return Sparse_Matrix;
   function Mult_M_SV (A : in Sparse_Matrix;
		       X : in Sparse_Vector) return Sparse_Vector
     with Pre => N_Col (A) = Length (X);
   function Permute_By_Col (Mat : in Sparse_Matrix;
			    P   : in Int_Array) return Sparse_Matrix;
   function Permute (Mat : in Sparse_Matrix;
		     P   : in Int_Array;
		     By  : in Permute_By_Type := Column) return Sparse_Matrix;
   procedure Transposed (Mat : in out Sparse_Matrix);
   

   function "-" (X : in Sparse_Matrix) return Sparse_Matrix;
   ---------- In Binary Form -----------------------------------------------
   function "*" (Left  : in Real;
		 Right : in Sparse_Matrix) return Sparse_Matrix;
   function "*" (Left  : in Sparse_Matrix;
		 Right : in Real) return Sparse_Matrix is (Right * Left);
   
   function "*" (Left, Right : in Sparse_Vector) return Sparse_Matrix;
   
   function "+" (Left, Right : in Sparse_Matrix) return Sparse_Matrix renames Plus2;
   function "-" (Left, Right : in Sparse_Matrix) return Sparse_Matrix renames Minus;
   function "*" (Left, Right : in Sparse_Matrix) return Sparse_Matrix renames Mult;
   function "*" (A : in Sparse_Matrix;
		 X : in Sparse_Vector) return Sparse_Vector renames Mult_M_SV;
   function "*" (A : in Sparse_Matrix;
		 X : in Real_Vector) return Sparse_Vector;
   
   function "and" (Left, Right : in Sparse_Matrix) return Sparse_Matrix renames Kronecker;
   function "or"  (Left, Right : in Sparse_Matrix) return Sparse_Matrix renames Direct_Sum;
   
   function "and" (Left  : in Sparse_Matrix;
		   Right : in Real_Matrix)   return Sparse_Matrix is (Left and Sparse (Right));
   function "and" (Left  : in Real_Matrix;
		   Right : in Sparse_Matrix) return Sparse_Matrix is (Sparse (Left) and Right);
   function "or" (Left  : in Sparse_Matrix;
		  Right : in Real_Matrix)    return Sparse_Matrix is (Left or Sparse (Right));
   function "or" (Left  : in Real_Matrix;
		  Right : in Sparse_Matrix)  return Sparse_Matrix is (Sparse (Left) or Right);

   function Remove_1stN (A : in Sparse_Matrix;
			 N : in Pos) return Sparse_Matrix;
     
   ------- File Readers ---------------------------------------------------
   function Read_Sparse_Triplet (File_Name : in String;
				 Offset	   : in Integer := 0) return Sparse_Matrix;
   
   
   
   procedure Testing_Stuff (A : in Sparse_Matrix);
private
   
   --  function BiCGSTAB (A   : in     Sparse_Matrix;
   --  		      B   : in     RVector;
   --  		      X0  : in     RVector;
   --  		      Err :    out Real;
   --  		      Tol : in     Real	    := 1.0e-10) return RVector;
   procedure Triplet_To_Matrix (Result :    out Sparse_Matrix;
				I      : in     IVector;
				J      : in     IVector;
				X      : in     RVector;
				N_Row  : in     Pos	      := 0;
				N_Col  : in     Pos	      := 0);
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   -------- Essential Tools -----------------------------------------
   procedure Cumulative_Sum (Item : in out Int_Array);
   procedure Remove_Duplicates (Mat : in out Sparse_Matrix);
   procedure Compress (Mat : in out Sparse_Matrix);
   -- Convert : goes from CSR to CSC or the reverse
   procedure Convert (Mat : in out Sparse_Matrix);
   
   ---- Define Matrix type -----------------------------------------
   --  type Sparse_Matrix is tagged
   --     record
   --  	 Format : Sparse_Matrix_Format := CSC;
   --  	 N_Row  : Pos := 0;
   --  	 N_Col  : Pos := 0;
   --  	 X      : RVector;
   --  	 I      : IVector;
   --  	 P      : IVector;
   --     end record;
   
   
   
   procedure Scatter (A	   : in     Sparse_Matrix;
		      J	   : in     Integer;
		      β	   : in     Real;
		      W	   : in out Int_Array;
		      X	   : in out Real_Vector;
		      Mark : in     Integer;
		      C	   : in out Sparse_Matrix;
		      Nz   : in out Integer);
   

   procedure To_Triplet (A     : in     Sparse_Matrix;
			 I     :    out IVector;
			 J     :    out IVector;
			 X     :    out RVector;
			 N_Row :    out Pos;
			 N_Col :    out Pos);
   
   type Sparse_Matrix is 
      record
	 Format : Sparse_Matrix_Format := CSC;
	 N_Row  : Pos := 0;
	 N_Col  : Pos := 0;
	 X      : RVector;
	 I      : IVector;
	 P      : IVector;
      end record;

end Numerics.Sparse_Matrices;
