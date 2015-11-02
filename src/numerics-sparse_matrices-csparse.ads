package Numerics.Sparse_Matrices.CSparse is
   
   ------- Define Matrix --------------------------------------------
   type LU_Type is private;
   
   ------ Basic Getter function ------------------------------------
   function N_Col (LU : in LU_Type) return Pos;
   
   ------ LU Decomposition -----------------------------------------
   function LU_Decomposition (Mat : in Sparse_Matrix;
			      Tol : in Real   := 1.0e-12) return LU_Type
     with Pre => Is_Valid (Mat) and Is_Square_Matrix (Mat);
   
   
   function Solve (LU : in LU_Type;
		   B  : in Real_Array) return Real_Array
     with Pre => N_Col (LU) = B'Length;
   function Solve (LU  : in LU_Type;
		   B   : in Sparse_Vector;
		   Tol : in Real	  := 1.0e-20) return Sparse_Vector
     with Pre => N_Col (LU) = Length (B);
   
   
   procedure Free (LU : in out LU_Type);
private
   ------- Define pointer packages ----------------------------------
   package Real_Ptrs is new C.Pointers (Index              => Nat,
					Element            => Real,
					Element_Array      => Real_Array,
					Default_Terminator => 0.0);
   package Int_Ptrs is new C.Pointers (Index              => Nat,
				       Element            => Int,
				       Element_Array      => Int_Array,
				       Default_Terminator => 0);

   ---- Define CS type -----------------------------------------
   type Sparse_Type is
      record
   	 Nzmax	: Pos	:= 0;
	 M	: Pos	:= 0;
	 N	: Pos	:= 0;
   	 I	: Int_Ptrs.Pointer;
   	 P	: Int_Ptrs.Pointer;
   	 X	: Real_Ptrs.Pointer;
	 Nz	: Pos	:= 0;
      end record with Convention => C;
   type Sparse_Ptr is access Sparse_Type with Convention => C;
   
   type Symbolic_Type is
      record
	 Pinv		: Int_Ptrs.Pointer;
	 Q		: Int_Ptrs.Pointer;
	 Parent		: Int_Ptrs.Pointer;
	 Cp		: Int_Ptrs.Pointer;
	 Leftmost	: Int_Ptrs.Pointer;
	 M2		: Int;
	 Lnz		: Real;
	 Unz		: Real;
      end record with Convention => C;
   type Symbolic_Ptr is access Symbolic_Type with Convention => C;
   
   type Numeric_Type is
      record
	 L	: Sparse_Ptr;
	 U	: Sparse_Ptr;
	 Pinv	: Int_Ptrs.Pointer;
	 B	: Real_Ptrs.Pointer;
      end record with Convention => C;
   type Numeric_Ptr  is access Numeric_Type with Convention => C;
   
   type LU_Type is
      record
	 Symbolic : Symbolic_Ptr;
	 Numeric  : Numeric_Ptr;
	 NCol     : Pos := 0;
      end record;
   
   ----------- Testing functions -------------------------------------------
   function Is_Valid (P	: in Real_Ptrs.Pointer;
		      N	: in Pos) return Boolean;

   
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
   
   function Solve_CS (N_Col : in     Int;
		      S	    : in     Symbolic_Ptr;
		      N	    : in     Numeric_Ptr;
		      B	    : in     Real_Array;
		      Err   :    out C.int) return Real_Ptrs.Pointer
     with Import => True, Convention => C, External_Name => "solve_cs";
   
   function Free (Sparse : in Sparse_Ptr) return Sparse_Ptr
     with Import => True, Convention => C, External_Name => "cs_dl_spfree";
   
   function Free (Symbolic : in Symbolic_Ptr) return Symbolic_Ptr
     with Import => True, Convention => C, External_Name => "cs_dl_sfree";
   
   function Free (Numeric : in Numeric_Ptr) return Numeric_Ptr
     with Import => True, Convention => C, External_Name => "cs_dl_nfree";
   
   procedure Print_Sparse (Sparse : in Sparse_Ptr)
     with Import => True, Convention => C, External_Name => "print_cs";
   
   function To_Sparse (Mat : in Sparse_Matrix) return Sparse_Ptr;
   
   function Solve (LU : in LU_Type;
		   B  : in Real_Array) return Real_Ptrs.Pointer
     with Post => Is_Valid (Solve'Result, B'Length);

   
   
end Numerics.Sparse_Matrices.CSparse;
