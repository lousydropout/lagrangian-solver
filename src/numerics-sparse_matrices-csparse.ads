with Interfaces.C, Interfaces.C.Pointers;
package Numerics.Sparse_Matrices.CSparse is
   
   ------- Define Matrix --------------------------------------------
   type LU_Type is private;
   
   ------ Basic Getter function ------------------------------------
   function N_Col (LU : in LU_Type) return Pos;
   
   ------ LU Decomposition -----------------------------------------
   function LU_Decomposition (Mat : in Sparse_Matrix;
			      Tol : in Real   := 1.0e-20) return LU_Type
     with Pre => Is_Valid (Mat) and Is_Square_Matrix (Mat);
   
   function Solve (A   : in Sparse_Matrix;
		   B   : in Sparse_Vector;
		   Tol : in Real	  := 1.0e-20) return Sparse_Vector;
   function Solve (LU  : in LU_Type;
		   B   : in Sparse_Vector;
		   Tol : in Real	  := 1.0e-20) return Sparse_Vector
     with Pre => N_Col (LU) = Length (B);
   
   procedure Free (LU : in out LU_Type);
   
private
   
   package C renames Interfaces.C;
   
   type Creal is new C.double range C.double'First .. C.double'Last;
   type Cint  is new C.long   range C.long'First   .. C.long'Last;
   subtype Cpos is Cint        range 0              .. Cint'Last;
   subtype Cnat is Cint        range 1              .. Cint'Last;
   
   type Creal_Array is array (Cnat range <>) of aliased Creal with Convention => C;
   type Cint_Array is array (Cnat range <>) of aliased Cint with Convention => C;
   
   ------- Define pointer packages ----------------------------------
   package Creal_Ptrs is new C.Pointers (Index              => Cnat,
					 Element            => Creal,
					 Element_Array      => Creal_Array,
					 Default_Terminator => 0.0);
   package Cint_Ptrs is new C.Pointers (Index              => Cnat,
					Element            => Cint,
					Element_Array      => Cint_Array,
					Default_Terminator => 0);

   ---- Define CS type -----------------------------------------
   type Sparse_Type is
      record
   	 Nzmax	: Cpos	:= 0;
	 M	: Cpos	:= 0;
	 N	: Cpos	:= 0;
   	 I	: Cint_Ptrs.Pointer;
   	 P	: Cint_Ptrs.Pointer;
   	 X	: Creal_Ptrs.Pointer;
	 Nz	: Cpos	:= 0;
      end record with Convention => C;
   type Sparse_Ptr is access Sparse_Type with Convention => C;
   
   type Symbolic_Type is
      record
	 Pinv		: Cint_Ptrs.Pointer;
	 Q		: Cint_Ptrs.Pointer;
	 Parent		: Cint_Ptrs.Pointer;
	 Cp		: Cint_Ptrs.Pointer;
	 Leftmost	: Cint_Ptrs.Pointer;
	 M2		: Cint;
	 Lnz		: Creal;
	 Unz		: Creal;
      end record with Convention => C;
   type Symbolic_Ptr is access Symbolic_Type with Convention => C;
   
   type Numeric_Type is
      record
	 L	: Sparse_Ptr;
	 U	: Sparse_Ptr;
	 Pinv	: Cint_Ptrs.Pointer;
	 B	: Creal_Ptrs.Pointer;
      end record with Convention => C;
   type Numeric_Ptr  is access Numeric_Type with Convention => C;
   
   type LU_Type is
      record
	 Symbolic : Symbolic_Ptr;
	 Numeric  : Numeric_Ptr;
	 NCol     : Cpos := 0;
      end record;
   
   ----------- Testing functions -------------------------------------------
   function Is_Valid (P	: in Creal_Ptrs.Pointer;
		      N	: in Cpos) return Boolean;

   
   ------------ C functions -------------------------------------------------
   function From_Arrays (M  : in Cint;
			 N  : in Cint;
			 Nz : in Cint;
			 I  : in Cint_Array;
			 J  : in Cint_Array;
			 X  : in Creal_Array) return Sparse_Ptr
     with Import => True, Convention => C, External_Name => "from_arrays";
   
   function To_CS (M	 : in Cint;
		   N	 : in Cint;
		   Nzmax : in Cint;
		   I	 : in Cint_Array;
		   P	 : in Cint_Array;
		   X	 : in Creal_Array) return Sparse_Ptr
     with Import => True, Convention => C, External_Name => "to_cs";
   
   function CS_Sqr (A	 : in Cint    := 0;
		    Prob : in Sparse_Ptr;
		    B	 : in Cint    := 0) return Symbolic_Ptr
     with Import => True, Convention => C, External_Name => "cs_dl_sqr";
   
   function CS_LU (Prob	: in Sparse_Ptr;
		   S	: in Symbolic_Ptr;
		   Tol	: in Creal    := 1.0e-15) return Numeric_Ptr
     with Import => True, Convention => C, External_Name => "cs_dl_lu";
   
   function Solve_CS (N_Col : in     Cint;
		      S	    : in     Symbolic_Ptr;
		      N	    : in     Numeric_Ptr;
		      B	    : in     Creal_Array;
		      Err   :    out C.int) return Creal_Ptrs.Pointer
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
		   B  : in Creal_Array) return Creal_Ptrs.Pointer
     with Post => Is_Valid (Solve'Result, B'Length);
   
   function Solve (LU : in LU_Type;
		   B  : in Creal_Array) return Creal_Array
     with Pre => N_Col (LU) = B'Length;
   function To_Array (X : in Real_Vector) return Creal_Array;
   function To_Array (X : in Int_Vector) return Cint_Array;   
   
end Numerics.Sparse_Matrices.CSparse;
