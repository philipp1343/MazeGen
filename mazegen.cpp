#include <cstdlib>
#include <functional>
#include <iostream>
#include <ostream>
#include <stack>
#include <set>
#include <vector>
#include <iosfwd>
#include <ctime>

class Cell {
public:
  Cell();

  bool hasWestWall() const;

  bool hasEastWall() const;

  bool hasNorthWall() const;

  bool hasSouthWall() const;

  bool isOnPath() const;

  void removeWestWall();

  void removeEastWall();

  void removeNorthWall();

  void removeSouthWall();

  void putOnPath();

  void removeFromPath();

private:
  bool m_hasWestWall;
  bool m_hasEastWall;
  bool m_hasNorthWall;
  bool m_hasSouthWall;
  bool m_isOnPath;
};

Cell::Cell()
  : m_hasWestWall(true)
  , m_hasEastWall(true)
  , m_hasNorthWall(true)
  , m_hasSouthWall(true)
  , m_isOnPath(false)
{
}

bool Cell::hasWestWall() const
{
  return m_hasWestWall;
}

bool Cell::hasEastWall() const
{
  return m_hasEastWall;
}

bool Cell::hasNorthWall() const
{
  return m_hasNorthWall;
}

bool Cell::hasSouthWall() const
{
  return m_hasSouthWall;
}

bool Cell::isOnPath() const
{
  return m_isOnPath;
}

void Cell::removeWestWall()
{
  m_hasWestWall = false;
}

void Cell::removeEastWall()
{
  m_hasEastWall = false;
}

void Cell::removeNorthWall()
{
  m_hasNorthWall = false;
}

void Cell::removeSouthWall()
{
  m_hasSouthWall = false;
}

void Cell::putOnPath()
{
  m_isOnPath = true;
}

void Cell::removeFromPath()
{
  m_isOnPath = false;
}
/*
 ^^^^^^^^^
 cell file
 */
class Position {
public:
  Position(int rowIndex, int columnIndex);

  int rowIndex() const;

  int columnIndex() const;

private:
  int m_rowIndex;
  int m_columnIndex;
};

Position::Position(int rowIndex, int columnIndex)
  : m_rowIndex(rowIndex), m_columnIndex(columnIndex)
{
}

int Position::rowIndex() const
{
  return m_rowIndex;
}

int Position::columnIndex() const
{
  return m_columnIndex;
}

bool operator==(Position lhs, Position rhs)
{
  return lhs.rowIndex() == rhs.rowIndex()
         && lhs.columnIndex() == rhs.columnIndex();
}

bool operator<(Position lhs, Position rhs)
{
  if (lhs.rowIndex() == rhs.rowIndex()) {
    return lhs.columnIndex() < rhs.columnIndex();
  }

  return lhs.rowIndex() < rhs.rowIndex();
}

/*
 ^^^^^^^^^^^^^
 position file
*/
void eraseIf(
  std::vector<Position>&               vector,
  const std::function<bool(Position)>& predicate)
{
  std::vector<Position> newVector;

  for (Position position : vector) {
    if (!predicate(position)) {
      newVector.push_back(position);
    }
  }

  vector = newVector;
}

/*
 ^^^^^^^^^^^^^
 erase if file
 */
namespace {
std::vector<Position> fetchNeighbors(Position position)
{
  std::vector<Position> neighbors;
  neighbors.push_back(Position(position.rowIndex() + 1, position.columnIndex()));
  neighbors.push_back(Position(position.rowIndex(), position.columnIndex() + 1));
  neighbors.push_back(Position(position.rowIndex(), position.columnIndex() - 1));
  neighbors.push_back(Position(position.rowIndex() - 1, position.columnIndex()));
  return neighbors;
}
} // anonymous namespace

std::vector<Position> fetchExistantNeighbors(
  Position    position,
  std::size_t rowCount,
  std::size_t columnCount)
{
  std::vector<Position> neighbors = fetchNeighbors(position);
  eraseIf(neighbors, [rowCount, columnCount](Position pos) {
    return pos.columnIndex() < 0 || pos.rowIndex() < 0
           || pos.columnIndex() >= columnCount || pos.rowIndex() >= rowCount;
  });
  return neighbors;
}
/*
 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 fetch existant neighbours file
 */

class Maze {
public:
  friend std::ostream& operator<<(std::ostream& os, const Maze& maze);

  Maze(int rows, int columns, int seed);

 
  bool isCellVisited(std::size_t rowIndex, std::size_t columnIndex) const;

  
  void visitCell(std::size_t rowIndex, std::size_t columnIndex);

  
  void unvisitCell(std::size_t rowIndex, std::size_t columnIndex);

  
  bool isPassageFree(Position from, Position to) const;

  std::size_t rowCount() const;

  std::size_t columnCount() const;

private:
  std::size_t                    m_rowCount;
  std::size_t                    m_columnCount;
  std::vector<std::vector<Cell>> m_matrix;
};


namespace {
std::vector<Position> fetchVisitableNeighbors(
  Position                            position,
  std::size_t                         rowCount,
  std::size_t                         columnCount,
  const std::set<Position>& visited)
{
  std::vector<Position> neighbors
    = fetchExistantNeighbors(position, rowCount, columnCount);
  eraseIf(
    neighbors, [&visited](Position pos) { return visited.count(pos) != 0; });
  return neighbors;
}

bool isNorthernNeighbor(Position source, Position target)
{
  return (target.rowIndex() == (source.rowIndex() - 1))
         && (source.columnIndex() == target.columnIndex());
}

bool isSouthernNeighbor(Position source, Position target)
{
  return (target.rowIndex() == (source.rowIndex() + 1))
         && (source.columnIndex() == target.columnIndex());
}

bool isWesternNeighbor(Position source, Position target)
{
  return (target.columnIndex() == (source.columnIndex() - 1))
         && (source.rowIndex() == target.rowIndex());
}

bool isEasternNeighbor(Position source, Position target)
{
  return (target.columnIndex() == (source.columnIndex() + 1))
         && (source.rowIndex() == target.rowIndex());
}

void breakDownWall(
  std::vector<std::vector<Cell>>& matrix,
  Position                        source,
  Position                        target)
{
  if (isNorthernNeighbor(source, target)) {
    matrix[source.rowIndex()][source.columnIndex()].removeNorthWall();
    matrix[target.rowIndex()][target.columnIndex()].removeSouthWall();
  }
  else if (isSouthernNeighbor(source, target)) {
    matrix[source.rowIndex()][source.columnIndex()].removeSouthWall();
    matrix[target.rowIndex()][target.columnIndex()].removeNorthWall();
  }
  else if (isWesternNeighbor(source, target)) {
    matrix[source.rowIndex()][source.columnIndex()].removeWestWall();
    matrix[target.rowIndex()][target.columnIndex()].removeEastWall();
  }
  else if (isEasternNeighbor(source, target)) {
    matrix[source.rowIndex()][source.columnIndex()].removeEastWall();
    matrix[target.rowIndex()][target.columnIndex()].removeWestWall();
  }
}

void recursiveBacktracker(std::vector<std::vector<Cell>>& matrix, int seed)
{
  std::srand(seed);
  const std::size_t    cellCount = matrix.size() * matrix.front().size();
  std::stack<Position> stack;
  std::set<Position> visited;
  stack.push(Position{0, 0});
  visited.insert(Position{0, 0});

  while (visited.size() != cellCount) {
    const Position              position  = stack.top();
    const std::vector<Position> neighbors = fetchVisitableNeighbors(
      position, matrix.size(), matrix.front().size(), visited);

    if (neighbors.empty()) {
      stack.pop();
      continue;
    }

    const int      index          = std::rand() % neighbors.size();
    const Position pickedNeighbor = neighbors[index];
    visited.insert(pickedNeighbor);
    stack.push(pickedNeighbor);
    breakDownWall(matrix, position, pickedNeighbor);
  }
}

void printNorthWall(std::ostream& os, const Cell& cell)
{
  os << '+';

  if (cell.hasNorthWall()) {
    os << "---";
  }
  else {
    os << "   ";
  }
}

void printSouthWall(std::ostream& os, const Cell& cell)
{
  os << '+';

  if (cell.hasSouthWall()) {
    os << "---";
  }
  else {
    os << "   ";
  }
}

void printWestWall(std::ostream& os, const Cell& cell)
{
  if (cell.hasWestWall()) {
    os << '|';
  }
  else {
    os << ' ';
  }
}

void printEastWall(std::ostream& os, const Cell& cell)
{
  if (cell.hasEastWall()) {
    os << '|';
  }
  else {
    os << ' ';
  }
}

void printMatrix(std::ostream& os, const std::vector<std::vector<Cell>>& matrix)
{
  const std::size_t rows       = matrix.size();
  const std::size_t columns    = matrix.front().size();
  const std::size_t lastRow    = rows - 1;
  const std::size_t lastColumn = columns - 1;

  for (std::size_t row = 0; row < rows; ++row) {
    
    for (std::size_t column = 0; column < columns; ++column) {
      const Cell& cell = matrix[row][column];
      printNorthWall(os, cell);

      if (column == lastColumn) {
        os << "+\n";
      }
    }

   
    for (std::size_t column = 0; column < columns; ++column) {
      const Cell& cell = matrix[row][column];
      printWestWall(os, cell);

      if (cell.isOnPath()) {
        os << " . ";
      }
      else {
        os << "   ";
      }

      if (column == lastColumn) {
        printEastWall(os, cell);
        os << '\n';
      }
    }

    
    if (row == lastRow) {
      for (std::size_t column = 0; column < columns; ++column) {
        const Cell& cell = matrix[row][column];
        printSouthWall(os, cell);

        if (column == lastColumn) {
          os << "+\n";
        }
      }
    }
  }
}
}

std::ostream& operator<<(std::ostream& os, const Maze& maze)
{
  const std::vector<std::vector<Cell>>& matrix = maze.m_matrix;
  printMatrix(os, matrix);
  return os;
}

Maze::Maze(int rows, int columns, int seed)
  : m_rowCount(rows)
  , m_columnCount(columns)
  , m_matrix(rows, std::vector<Cell>(columns))
{
  recursiveBacktracker(m_matrix, seed);
}

bool Maze::isCellVisited(std::size_t rowIndex, std::size_t columnIndex) const
{
  return m_matrix[rowIndex][columnIndex].isOnPath();
}

void Maze::visitCell(std::size_t rowIndex, std::size_t columnIndex)
{
  m_matrix[rowIndex][columnIndex].putOnPath();
}

void Maze::unvisitCell(std::size_t rowIndex, std::size_t columnIndex)
{
  m_matrix[rowIndex][columnIndex].removeFromPath();
}

bool Maze::isPassageFree(Position from, Position to) const
{
  if (isNorthernNeighbor(from, to)) {
    return !m_matrix[from.rowIndex()][from.columnIndex()].hasNorthWall();
  }

  if (isSouthernNeighbor(from, to)) {
    return !m_matrix[from.rowIndex()][from.columnIndex()].hasSouthWall();
  }

  if (isWesternNeighbor(from, to)) {
    return !m_matrix[from.rowIndex()][from.columnIndex()].hasWestWall();
  }

  if (isEasternNeighbor(from, to)) {
    return !m_matrix[from.rowIndex()][from.columnIndex()].hasEastWall();
  }

  return false;
}

std::size_t Maze::rowCount() const
{
  return m_rowCount;
}

std::size_t Maze::columnCount() const
{
  return m_columnCount;
}
/*
 ^^^^^^^^^^^^^
 maze file
 */

namespace {
std::vector<Position> directNeighborsThatCanBeReached(
  Position    from,
  const Maze& maze)
{
  std::vector<Position> neighbors
    = fetchExistantNeighbors(from, maze.rowCount(), maze.columnCount());
  eraseIf(neighbors, [from, &maze](Position to) {
    return !maze.isPassageFree(from, to);
  });
  return neighbors;
}
} // anonymous namespace

bool findPath(Maze& maze, Position from, Position to)
{
  maze.visitCell(from.rowIndex(), from.columnIndex());

  if (from == to) {
    return true;
  }

  std::vector<Position> neighbors = directNeighborsThatCanBeReached(from, maze);

  for (Position n : neighbors) {
    if (!maze.isCellVisited(n.rowIndex(), n.columnIndex())) {
      if (findPath(maze, n, to)) {
        return true;
      }
    }
  }

  maze.unvisitCell(from.rowIndex(), from.columnIndex());
  return false;
}

/*
 ^^^^^^^^^^^^^^
 find path file
 */




int createSeed(int argc, char** argv)
{
  const int seedIndex = 3;

  if (seedIndex < argc) {
    const char* seedString = argv[seedIndex];
    const int   seed       = std::atoi(seedString);
    return seed;
  }

  return static_cast<int>(std::time(nullptr));
}
/*
 ^^^^^^^^^^^^^^^^
 create seed file
 */



int main(int argc, char* argv[])
{
  const int minimumArgumentCount = 3;
  const int maximumArgumentCount = 4;

  if ((argc < minimumArgumentCount) || (argc > maximumArgumentCount)) {
    std::cerr << "Invalid argument count (" << (argc - 1) << ").\n";
    return EXIT_FAILURE;
  }

  const int      rowsIndex     = 1;
  const int      columnsIndex  = 2;
  const char*    rowsString    = argv[rowsIndex];
  const char*    columnsString = argv[columnsIndex];
  const int      rows          = std::atoi(rowsString);
  const int      columns       = std::atoi(columnsString);
  const int      seed          = createSeed(argc, argv);
  const Position from          = {0, 0};
  const Position to            = {rows - 1, columns - 1};

  Maze       maze(rows, columns, seed);
   findPath(maze, from, to);
  std::cout <<maze << '\n';
}

