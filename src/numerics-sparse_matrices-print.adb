separate (Numerics.Sparse_Matrices)

procedure Print (Mat : in Sparse_Matrix) is
   
   use Ada.Text_IO, Sparse_Matrix_Format_IO, Int_IO, Real_IO;
   
begin
   
   Put ("  (");
   Put (Mat.N_Row, Width => 0); Put (" x ");
   Put (Mat.N_Col, Width => 0);
   Put (") matrix in "); Put (Mat.Format); Put (" format.");
   New_Line;
   
   
   Put ("I: ");
   for I of Mat.I loop
      Put (", "); Put (I, Width => 4);
   end loop;
   New_Line; 
   
   Put ("P: ");
   for P of Mat.P loop
      Put (", "); Put (P, Width => 4);
   end loop;
   New_Line; 
   
   Put ("X: ");
   for X of Mat.X loop 
      Put (", "); Put (X, Aft => 3, Exp => 2, Fore => 3);
   end loop;
   New_Line;
   
end Print;
