1. First, the Point class is implemented to represent a point on the GIS system, and overloads its commonly used operator functions to facilitate subsequent use.

2. Implement the Quad tree:
(1) In order to realize this Quad tree, we must first implement its Node Node class. This Node class mainly contains the lower left, upper right and central nodes of this Node, while maintaining half the length and half width of the rectangle. At the same time, this node provides an interface function to judge whether it contains or intersects with other nodes.  (2)  Implement the QuadTree class, which is the concrete implementation of the Quad tree, which contains all the Point and the four sub-Quad trees.   The function of inserting   and querying a node is provided here, which is mainly realized by recursion. At the same time, this class also provides a clear interface, which is used to clear the state of the current class and prepare to build the Quad tree next time. The function of inserting and querying a node is provided here, which is mainly realized by recursion. At the same time, this class also provides a clear interface, which is used to clear the state of the current class and prepare to build the Quad tree next time.

4. Implement the Hash table:  Here, the probing is used to deal with hash conflicts. If it is found that the current hash position is occupied, the square is used to calculate the next hash position. At the same time, when the data capacity exceeds 70% size, the hash table will perform double size operation and rehash at the same time. This hash table provides a clear interface to clear the current hash table and prepare for the next hash operation.

5. Then implement the DMS and DEC classes, which use Strings and Floats to represent latitude and longitude information, and provide conversion functions between the two.

6. Implement GeoDataRecord, which represents a row record in the database;
The properties of this class correspond to the column information of the database one by one. At the same time, tolongstring function is provided to standardize long print.

7. Implement Cacheitem, which stores the time information of the accessed data item to prepare for maintaining the buffer pool

8. Implement Geodatabase, which implements database functions, using the Quadtree and HashMap implemented previously to maintain database information. It can read and write data. At the same time, this class implements the function of searching according to name and coordinate. This class also maintains the buffer pool and updates the cache item. In addition, point to Geofeature conversion is provided.

9. Implement commands and tokenizer class aid classes to parse commands

10. Finally, the GIS class is implemented, which accepts three files: database file, script file and log file. It processes the instructions in script file one by one, and writes the data and log in databse file and log file.
