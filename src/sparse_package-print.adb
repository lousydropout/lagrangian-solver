separate (Sparse_Package)

procedure Print (Mat : in Matrix) is
   
   use Ada.Text_IO, Matrix_Format_IO, Int_IO, Real_IO;
   
begin
   
   Put ("  (");
   Put (Mat.N_Row, Width => 0); Put (" x ");
   Put (Mat.N_Col, Width => 0);
   Put (") matrix in "); Put (Mat.Format); Put (" format.");
   New_Line;
   
   
   Put ("I: ");
   for I of Mat.I loop
      Put (", "); Put (I, Width => 6);
   end loop;
   New_Line; 
   
   Put ("P: ");
   for P of Mat.P loop
      Put (", "); Put (P, Width => 6);
   end loop;
   New_Line; 
   
   Put ("X: ");
   for X of Mat.X loop 
      Put (", "); Put (X, Aft => 1, Exp => 2, Fore => 3);
   end loop;
   New_Line;
   
end Print;
