with Ada.Text_IO; use Ada.Text_IO;
with Sparse_Package; use Sparse_Package; 

procedure Sparse_Test is
   N : Int := 2;
   
   I1 : Int_Array  := (3,   2,   1);
   J1 : Int_Array  := (1,   3,   2);
   X1 : Real_Array := (1.234, 2.345, 2.789);
   Left : Matrix := Triplet_To_Matrix (I1, J1, X1, 3, 3);
   
   
   I2 : Int_Array  := (1,   2, 2);
   J2 : Int_Array  := (1,   2, 1);
   X2 : Real_Array := (1.5, 3.0, 1.0);
   Right : Matrix 
     := Triplet_To_Matrix (I => I2, 
			   J => J2, 
			   X => X2, 
			   N_Row => N, 
			   N_Col => N + 1);
   
   I3 : Int_Array  := (5, 6, 6, 1, 2, 2, 3, 4, 4);
   J3 : Int_Array  := (1, 1, 2, 4, 4, 5, 7, 7, 8);
   X3 : Real_Array := (1.5, 1.0, 3.0, 3.0, 2.0, 6.0, 3.0, 2.0, 6.0);
   Right2 : Matrix := Triplet_To_Matrix (I3, J3, X3, 6, 9);
   
   X : Real_Array (1 .. N) := (others => 0.0);
   Vec, X0 : Real_Vector;
   
   Result : Sparse_Ptr;
   
   use Real_Ptrs;
   function Create (I : in Int) return Pointer
     with Import => True, Convention => C, External_Name => "dvec";
   P : Pointer := Create (5);
   
   
   
   --  Ar : Real_Array := Value (P, 5);
   Ar : array (1 .. 5) of aliased Real
     with Convention => C, Address => P.all'Address;
begin
   Put_Line (Real'Image (P.all));
   for K in 1 .. 5 loop
      Put_Line (Real'Image (Ar (K)));
      --  Put_Line (Real'Image (Value (P, 5) (K)));
   end loop;
   New_Line;
   
   Left.Print;
   New_Line;
   Put_Line ("Number of Elements = " & Int'Image (Number_Of_Elements (Left)));
   Result := To_Sparse (Left);
   Print_Sparse (Result);
end Sparse_Test;
