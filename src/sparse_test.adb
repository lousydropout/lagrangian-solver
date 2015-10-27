with Ada.Text_IO; use Ada.Text_IO;
with Sparse_Package; use Sparse_Package; 

procedure Sparse_Test is
   use Real_Ptrs, Real_IO;
   function Create (I : in Int) return Pointer
     with Import => True, Convention => C, External_Name => "dvec";
   P : Pointer := Create (5);
   Ar : array (1 .. 5) of aliased Real
     with Convention => C, Address => P.all'Address;
   
   
   
   
   N : Int := 2;
   
   I1 : Int_Array  := (3,   2,   1);
   J1 : Int_Array  := (1,   3,   2);
   X1 : Real_Array := (1.234, 2.345, 2.789);
   Left : Matrix := Triplet_To_Matrix (I1, J1, X1, 3, 3);
   
   X : Real_Array (1 .. N) := (others => 0.0);
   Vec, X0 : Real_Vector;
   
   Result : Sparse_Ptr;
   
begin
   Put_Line (Real'Image (P.all));
   for K of Ar loop
      Put (K); New_Line;
   end loop;
   New_Line;
   
   Left.Print;
   New_Line;
   Put_Line ("Number of Elements = " & Int'Image (Number_Of_Elements (Left)));
   Result := To_Sparse (Left);
   Print_Sparse (Result);
end Sparse_Test;
