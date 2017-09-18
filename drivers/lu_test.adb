with Numerics, Numerics.Dense_Matrices, Ada.Text_IO, Auto_Differentiation_Dense;
use  Numerics, Numerics.Dense_Matrices, Ada.Text_IO;

procedure LU_Test is
   use Real_IO, Int_IO;
   N : constant Nat := 5;
   package AD is new Auto_Differentiation_Dense (N); use AD;
   
   Example_1 : constant Real_Matrix := ((1.0, 3.0, 5.0),
					(2.0, 4.0, 7.0),
					(1.0, 1.0, 0.0));
   Inv, L1, U1 : Real_Matrix (Example_1'Range (1), Example_1'Range (2));
   P1     : Int_Array  (Example_1'Range (1));
   
   Example_2 : constant Real_Matrix := ((11.0, 9.0, 24.0, 2.0),
					(1.0, 5.0, 2.0, 6.0),
					(3.0, 17.0, 18.0, 1.0),
					(2.0, 5.0, 7.0, 1.0));
   L2, U2 : Real_Matrix (Example_2'Range (1), Example_2'Range (2));
   P2     : Int_Array  (Example_2'Range (1));
   
   Num : Integer;
   X : Real_Vector := (-1.0, 3.0, 2.0);
begin

   LU_Decomposition (A => Example_1,
		     P => P1,
		     L => L1,
		     U => U1);
   LU_Decomposition (A => Example_2,
		     P => P2,
		     L => L2,
		     U => U2);
   

   Put_Line ("Example 1: ");
   Print (Example_1);
   Put_Line ("P: ");
   for P of P1 loop
      Put (P, 3); Put (",  ");
   end loop;
   New_Line;
   Put_Line ("L: ");
   Print (L1);
   Put_Line ("U: ");
   Print (U1);
   Put ("Det = "); Put (Determinant (P1, L1, U1), Exp => 0); New_Line;
   Inv := Inverse (Example_1);
   Put_Line ("Inv: ");
   Print (Inv);
   Inv := Inv * Example_1;
   New_Line;
   Print (Inv);
   New_Line;
   X := Solve (P1, L1, U1, X);
   Put_Line ("X: ");
   for Item of X loop
      Put (Item, Aft => 3, Exp => 0); New_Line;
   end loop;
   New_Line;
   
   --  Put_Line ("Example 2: ");
   --  Print (Example_2);
   --  Put_Line ("P: ");
   --  for P of P2 loop
   --     Put (P, 3); Put (",  ");
   --  end loop;
   --  New_Line;
   --  Put_Line ("L: ");
   --  Print (L2);
   --  Put_Line ("U: ");
   --  Print (U2);
   --  Put ("Det = "); Put (Determinant (P2, L2, U2), Exp => 0); New_Line;
   --  Put ("Det = "); Put (Determinant (Example_2), Exp => 0); New_Line;
   --  New_Line;
   
   
   --  Put_Line ("----------------------------------------");
   --  Num := Number_Of_Swaps ((1, 2, 3, 4, 5));
   --  Put ("Number of swaps = "); Put (Num); New_Line;
   --  Put_Line ("----------------------------------------");
   
end LU_Test;
