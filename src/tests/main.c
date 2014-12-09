void connexityTests();
void convexSolverTests();
void dijkstraTests();
void kmeansClusteringTests();
void neighborsTests();

// Main function calling all unit tests
// If nothing gets printed, everything's OK
int main()
{
	test_connexity1();
	test_connexity2();

	test_convexSolver1();
	test_convexSolver2();

	test_dijkstra1();
	test_dijkstra2();

	test_kmeansClustering1();
	test_kmeansClustering2();

	test_neighbors1();
	test_neighbors2();

	return 0;
}
