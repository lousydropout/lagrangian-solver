with Ada.Containers.Vectors, Interfaces.C, Interfaces.C.Pointers;

package Sparse_Package is
   package C renames Interfaces.C;
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   -------- Define Enumeration types --------------------------------
   type Permute_By_Type is (Row, Column);
   type Matrix_Format   is (CSR, CSC, Triplet);
   -------- Define types (Real, Int, Pos, Nat) ----------------------
   type Real is new C.double range C.double'First .. C.double'Last;
   subtype Int is Integer range Integer'First .. Integer'Last;
   subtype Pos is Int     range 0 .. Int'Last;
   subtype Nat is Pos     range 1 .. Pos'Last;
   -------- Define array types
   type Real_Array is array (Nat range <>) of Real;
   type Int_Array  is array (Nat range <>) of Int;
   
   
   ------- Define Real_Vector and Int_Vector ------------------------
   package IV_Package is new Ada.Containers.Vectors (Nat, Int, "=");
   package RV_Package is new Ada.Containers.Vectors (Nat, Real, "=");
   subtype Int_Vector  is IV_Package.Vector;
   subtype Real_Vector is RV_Package.Vector;
   ------- Define Matrix --------------------------------------------
   type Matrix is tagged private;
   
   --- Print procedure
   procedure Print (Mat : in Matrix);  -- The only public procedure
   
   ------- Basic Getter Functions -----------------------------------
   function N_Row (Mat : in Matrix) return Pos;
   function N_Col (Mat : in Matrix) return Pos;
   function Max (Item : in Int_Array) return Int with Inline => True;
   function Abs_Max (Item : in Int_Array) return Int with Inline => True;
   function Max (Item : in Real_Array) return Real with Inline => True;
   function Abs_Max (Item : in Real_Array) return Real with Inline => True;

   
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
private
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   -------- Essential Tools -----------------------------------------
   function  Cumulative_Sum (Item : in Int_Array) return Int_Array;
   --  procedure Transposed (Mat : in out Matrix) with Inline => True;
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
      
end Sparse_Package;
