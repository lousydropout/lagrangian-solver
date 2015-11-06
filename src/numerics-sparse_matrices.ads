package Numerics.Sparse_Matrices is
   
   
   -------- Define Enumeration types --------------------------------
   type Permute_By_Type is (Row, Column);
   type Sparse_Matrix_Format   is (CSR, CSC, Triplet);
   package Sparse_Matrix_Format_IO is new Ada.Text_IO.Enumeration_IO (Sparse_Matrix_Format);
   
   ------- Define Matrix --------------------------------------------
   type Sparse_Matrix is tagged private;
   
   --- Print procedure ----------------------------------------------
   procedure Print (Mat : in Sparse_Matrix); 
   
   ------- Basic Getter Functions -----------------------------------
   function Norm2 (Item : in Sparse_Matrix) return Real;
   function N_Row (Mat : in Sparse_Matrix)  return Pos;
   function N_Col (Mat : in Sparse_Matrix)  return Pos;
   function Number_Of_Elements (X : in Sparse_Matrix) return Int;
   
   ------- Functions for Creating Sparse Matrices -------------------
   procedure Set_Diag (X  : in out Sparse_Matrix;
   		       To : in     Sparse_Vector)
     with Pre => Is_Square_Matrix (X) and X.N_Col = Length (To);
   function Diag (X : in Sparse_Matrix) return Sparse_Vector
     with Pre => Is_Square_Matrix (X);
   function Diag (X : in Sparse_Vector) return Sparse_Matrix;
   function Sparse (X : in Real_Matrix) return Sparse_Matrix;
   function Triplet_To_Matrix (I      : in Int_Array;
			       J      : in Int_Array;
			       X      : in Real_Array;
			       N_Row  : in Pos := 0;
			       N_Col  : in Pos := 0;
			       Format : in Sparse_Matrix_Format := CSC) return Sparse_Matrix
     with Pre => I'Length = J'Length and I'Length = X'Length;
   function Convert (Mat : in Sparse_Matrix) return Sparse_Matrix;
   procedure Add (Mat  : in out Sparse_Matrix;
		  I, J : in     Nat;
		  X    : in     Real)
     with Pre => I <= Mat.N_Row and J <= Mat.N_Col;
   procedure Set (Mat  : in out Sparse_Matrix;
		  I, J : in     Nat;
		  X    : in     Real)
     with Pre => I <= Mat.N_Row and J <= Mat.N_Col;
   
   
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
   function Transpose (Mat : in Sparse_Matrix) return Sparse_Matrix;
   function Plus (Left  : in Sparse_Matrix;
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
   

   ---------- In Binary Form -----------------------------------------------
   function "+" (Left, Right : in Sparse_Matrix) return Sparse_Matrix renames Plus;
   function "-" (Left, Right : in Sparse_Matrix) return Sparse_Matrix renames Minus;
   function "*" (Left, Right : in Sparse_Matrix) return Sparse_Matrix renames Mult;
   function "*" (A : in Sparse_Matrix;
		 X : in Sparse_Vector) return Sparse_Vector renames Mult_M_SV;
   
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

   
   ------- File Readers ---------------------------------------------------
   function Read_Sparse_Triplet (File_Name : in String;
				 Offset	   : in Int    := 0) return Sparse_Matrix;
   
   
   
private
   
   --  function BiCGSTAB (A   : in     Sparse_Matrix;
   --  		      B   : in     Real_Vector;
   --  		      X0  : in     Real_Vector;
   --  		      Err :    out Real;
   --  		      Tol : in     Real	    := 1.0e-10) return Real_Vector;
   procedure Triplet_To_Matrix (Result :    out Sparse_Matrix;
				I      : in     Int_Vector;
				J      : in     Int_Vector;
				X      : in     Real_Vector;
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
   type Sparse_Matrix is tagged
      record
	 Format : Sparse_Matrix_Format;
	 N_Row  : Pos := 0;
	 N_Col  : Pos := 0;
	 X      : Real_Vector;
	 I      : Int_Vector;
	 P      : Int_Vector;
      end record;
   
   
   
   procedure Scatter (A	   : in     Sparse_Matrix;
		      J	   : in     Int;
		      β	   : in     Real;
		      W	   : in out Int_Array;
		      X	   : in out Real_Array;
		      Mark : in     Int;
		      C	   : in out Sparse_Matrix;
		      Nz   : in out Int);
   

end Numerics.Sparse_Matrices;
