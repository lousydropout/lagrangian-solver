with Ada.Containers.Vectors, Interfaces.C, Interfaces.C.Pointers, Ada.Text_IO;

package Sparse_Package is
   package C renames Interfaces.C;
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   -------- Define Enumeration types --------------------------------
   type Permute_By_Type is (Row, Column);
   type Matrix_Format   is (CSR, CSC, Triplet);
   
   -------- Define types (Real, Int, Pos, Nat) ----------------------
   type Real is new C.double range C.double'First .. C.double'Last;
   type Int  is new C.long range C.long'First .. C.long'Last;
   subtype Pos is Int     range 0 .. Int'Last;
   subtype Nat is Pos     range 1 .. Pos'Last;
   
   -------- Define array types
   type Real_Array is array (Nat range <>) of aliased Real with Convention => C;
   type Int_Array  is array (Nat range <>) of aliased Int  with Convention => C;
   
   ------- Define Real_IO and Int_IO packages ------------------------
   package Int_IO is new Ada.Text_IO.Integer_IO (Int);
   package Real_IO is new Ada.Text_IO.Float_IO (Real);
   package Matrix_Format_IO is new Ada.Text_IO.Enumeration_IO (Matrix_Format);
   ------- Define Real_Vector and Int_Vector packages ------------------------
   package IV_Package is new Ada.Containers.Vectors (Nat, Int, "=");
   package RV_Package is new Ada.Containers.Vectors (Nat, Real, "=");
   subtype Int_Vector  is IV_Package.Vector;
   subtype Real_Vector is RV_Package.Vector;
   
   ------- Define pointer packages ----------------------------------
   package Real_Ptrs is new C.Pointers (Index              => Nat,
					Element            => Real,
					Element_Array      => Real_Array,
					Default_Terminator => 0.0);
   package Int_Ptrs is new C.Pointers (Index              => Nat,
				       Element            => Int,
				       Element_Array      => Int_Array,
				       Default_Terminator => 0);
   
   ------- Define Matrix --------------------------------------------
   type Matrix        is tagged private;
   type LU_Type       is private;
   type Sparse_Type   is private; type Sparse_Ptr   is private;
   --- Print procedure
   procedure Print (Mat : in Matrix);  -- The only public procedure
   
   ------- Basic Getter Functions -----------------------------------
   function Norm2 (Item : in Matrix) return Real;
   function Norm2_RV (X : in Real_Vector) return Real;
   function Norm_RV (X : in Real_Vector) return Real;
   function N_Row (Mat : in Matrix) return Pos with Inline => True;
   function N_Col (Mat : in Matrix) return Pos with Inline => True;
   function Max_Int_Array (Item : in Int_Array) return Int with Inline => True;
   function Abs_Max_IA (Item : in Int_Array) return Int with Inline => True;
   function Max_Real_Array (Item : in Real_Array) return Real with Inline => True;
   function Abs_Max_RA (Item : in Real_Array) return Real with Inline => True;

   function Max (Item : in Int_Array) return Int renames Max_Int_Array;
   function Max (Item : in Real_Array) return Real renames Max_Real_Array;
   function Abs_Max (Item : in Int_Array) return Int renames Abs_Max_IA;
   function Abs_Max (Item : in Real_Array) return Real renames Abs_Max_RA;
   function Norm2 (X : in Real_Vector) return Real renames Norm2_RV;
   function Norm (X : in Real_Vector) return Real renames Norm_RV;
   
   function Length (X : in Real_Vector) return Int;
   function Number_Of_Elements (X : in Matrix) return Int;
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Functions for Creating Sparse Matrices -------------------
   function Triplet_To_Matrix (I      : in Int_Array;
			       J      : in Int_Array;
			       X      : in Real_Array;
			       N_Row  : in Pos := 0;
			       N_Col  : in Pos := 0;
			       Format : in Matrix_Format := CSC) return Matrix
     with Pre => I'Length = J'Length and I'Length = X'Length;
   
   function Convert (Mat : in Matrix) return Matrix;
   function Vectorize (I : in Int_Array;
		       X : in Real_Array) return Matrix
     with Pre => I'Length = X'Length;
   
   
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Testing Functions -----------------------------------
   function Is_Col_Vector (A : in Matrix) return Boolean;
   function Is_Square_Matrix (A : in Matrix) return Boolean;
   function Has_Same_Dimensions (Left, Right : in Matrix) return Boolean;
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Matrix operations -----------------------------------
   function Eye (N : in Nat) return Matrix;
   function Zero_Vector (N : in Nat) return Matrix;
   function Dot_Product (Left_I, Right_J : in Int_Array;
			 Left_X, Right_Y : in Real_Array) return Real;
   function Dot_Product_RV (X, Y : in Real_Vector) return Real;
   function Dot_Product (X, Y : in Real_Vector) return Real 
     renames Dot_Product_RV;
   function Transpose (Mat : in Matrix) return Matrix;
   function Plus (Left  : in Matrix;
		  Right : in Matrix) return Matrix
     with Pre => Has_Same_Dimensions (Left, Right);
   function Minus (Left  : in Matrix;
		   Right : in Matrix) return Matrix
     with Pre => Has_Same_Dimensions (Left, Right);
   function Mult (Left, Right : in Matrix) return Matrix
     with Pre => N_Col (Left) = N_Row (Right);
   function Mult_Int_Array (Left, Right : in Int_Array) return Boolean;
   function Kronecker (Left, Right : in Matrix) return Matrix;
   function Direct_Sum (Left, Right : in Matrix) return Matrix;
   function Mult_R_RV (Left  : in Real;
		       Right : in Real_Vector) return Real_Vector;
   function Mult_M_RV (Left  : in Matrix;
		       Right : in Real_Vector) return Real_Vector
     with Pre => N_Row (Left) = Pos (Right.Length);
   function Add_RV_RV (Left, Right : in Real_Vector) return Real_Vector;
   function Minus_RV_RV (Left, Right : in Real_Vector) return Real_Vector;
   function Permute_By_Col (Mat : in Matrix;
			    P   : in Int_Array) return Matrix;
   function Permute (Mat : in Matrix;
		     P   : in Int_Array;
		     By  : in Permute_By_Type := Column) return Matrix;
   procedure Transposed (Mat : in out Matrix);
   

   ---------- In Binary Form -----------------------------------------------
   function "+" (Left, Right : in Matrix) return Matrix renames Plus;
   function "-" (Left, Right : in Matrix) return Matrix renames Minus;
   function "*" (Left, Right : in Matrix) return Matrix renames Mult;
   function "*" (Left, Right : in Int_Array) return Boolean renames Mult_Int_Array;
   function "and" (Left, Right : in Matrix) return Matrix renames Kronecker;
   function "or" (Left, Right : in Matrix) return Matrix renames Direct_Sum;
   function "*" (Left  : in Real;
		 Right : in Real_Vector) return Real_Vector renames Mult_R_RV;
   function "*" (Left  : in Matrix;
		 Right : in Real_Vector) return Real_Vector renames Mult_M_RV;
   function "-" (Left, Right : in Real_Vector) return Real_Vector renames Minus_RV_RV;
   function "+" (Left, Right : in Real_Vector) return Real_Vector renames Add_RV_RV;
   
   
   
   
   
   function LU_Decomposition (Mat : in Matrix;
			      Tol : in Real   := 1.0e-12) return LU_Type;
   function Solve (LU : in LU_Type;
		   B  : in Real_Array) return Real_Ptrs.Pointer;
   function Solve (LU : in LU_Type;
		   B  : in Real_Array) return Real_Array;
   
   ----------------- Ada wrappers of C functions -------------------------------
   function To_Sparse (Mat : in Matrix) return Sparse_Ptr;
   procedure Print_Sparse (Sparse : in Sparse_Ptr)
     with Import => True, Convention => C, External_Name => "print_cs";
   
private
   
   function BiCGSTAB (A   : in     Matrix;
		      B   : in     Real_Vector;
		      X0  : in     Real_Vector;
		      Err :    out Real;
		      Tol : in     Real	    := 1.0e-10) return Real_Vector;
   

   ------------------------------------------------------------------
   ------------------------------------------------------------------
   -------- Essential Tools -----------------------------------------
   function  Cumulative_Sum (Item : in Int_Array) return Int_Array;
   procedure Remove_Duplicates (Mat : in out Matrix);
   procedure Compress (Mat : in out Matrix);
   -- Convert : goes from CSR to CSC or the reverse
   procedure Convert (Mat : in out Matrix);
   
   -- Vectorize & To_Array are needed in Triplet_To_Matrix
   function Vectorize (Item : in Real_Array) return Real_Vector;
   function Vectorize (Item : in Int_Array)  return Int_Vector;
   function To_Array (Item : in Real_Vector) return Real_Array;
   function To_Array (Item : in Int_Vector) return Int_Array;
   

   
   
   
   ---- Define Matrix type -----------------------------------------
   type Matrix is tagged
      record
	 Format : Matrix_Format;
	 N_Row  : Nat;
	 N_Col  : Nat;
	 X      : Real_Vector;
	 I      : Int_Vector;
	 P      : Int_Vector;
      end record;
   
   ---- Define CS type -----------------------------------------
   
   type Sparse_Type is
      record
   	 Nzmax, M, N : Pos;
   	 P, I : Int_Ptrs.Pointer;
   	 X : Real_Ptrs.Pointer;
	 Nz : Pos;
      end record with Convention => C;
   type Symbolic_Type is
      record
	 Pinv, Q, Parent, Cp, Leftmost : Int_Ptrs.Pointer;
	 M2 : Int;
	 Lnz, Unz : Real;
      end record with Convention => C;
   type Numeric_Type is
      record
	 L, U : Sparse_Ptr;
	 Pinv : Int_Ptrs.Pointer;
	 B : Real_Ptrs.Pointer;
      end record with Convention => C;
   type Symbolic_Ptr is access Symbolic_Type with Convention => C;
   type Numeric_Ptr  is access Numeric_Type  with Convention => C;
   
   type LU_Type is
      record
	 Symbolic : Symbolic_Ptr;
	 Numeric  : Numeric_Ptr;
	 NCol     : Pos;
      end record;
   
   type Sparse_Ptr   is access Sparse_Type   with Convention => C;
   
   
   
   ------------ C functions -------------------------------------------------
   function From_Arrays (M  : in Int;
			 N  : in Int;
			 Nz : in Int;
			 I  : in Int_Array;
			 J  : in Int_Array;
			 X  : in Real_Array) return Sparse_Ptr
     with Import => True, Convention => C, External_Name => "from_arrays";
   
   function To_CS (M	 : in Int;
		   N	 : in Int;
		   Nzmax : in Int;
		   I	 : in Int_Array;
		   P	 : in Int_Array;
		   X	 : in Real_Array) return Sparse_Ptr
     with Import => True, Convention => C, External_Name => "to_cs";
   
   function CS_Sqr (A	 : in Int    := 0;
		    Prob : in Sparse_Ptr;
		    B	 : in Int    := 0) return Symbolic_Ptr
     with Import => True, Convention => C, External_Name => "cs_dl_sqr";
   function CS_LU (Prob	: in Sparse_Ptr;
		   S	: in Symbolic_Ptr;
		   Tol	: in Real    := 1.0e-15) return Numeric_Ptr
     with Import => True, Convention => C, External_Name => "cs_dl_lu";
   function Solve_CS (N_Col : in Int;
   		      S : in Symbolic_Ptr;
   		      N : in Numeric_Ptr;
   		      B : in Real_Array) return Real_Ptrs.Pointer
     with Import => True, Convention => C, External_Name => "solve_cs";
   
   function Free (Sparse : in Sparse_Ptr) return Sparse_Ptr
      with Import => True, Convention => C, External_Name => "cs_dl_spfree";
   function Free (Symbolic : in Symbolic_Ptr) return Symbolic_Ptr
      with Import => True, Convention => C, External_Name => "cs_dl_sfree";
   function Free (Numeric : in Numeric_Ptr) return Numeric_Ptr
      with Import => True, Convention => C, External_Name => "cs_dl_nfree";
end Sparse_Package;
