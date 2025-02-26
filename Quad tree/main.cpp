#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <chrono>
#include <functional>
#include <list>
#include <set>
#include <unordered_set>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <unordered_map>

using namespace std;

vector<string> Split(const string &input, const string &delimiter, const bool skip_empty = false) {
    vector<string> output;
    auto first = input.begin();
    while (first != input.end()) {
        const auto second = find_first_of(first, input.end(), delimiter.begin(), delimiter.end());
        if (first != second || !skip_empty) {
            output.emplace_back(first, second);
        }
        if (second == input.end()) {
            break;
        }
        first = next(second);
    }
    return output;
}

// A Point class
class Point {
public:
    float x;
    float y;
    vector<size_t> indices;

public:
    explicit Point(const float x = 0, const float y = 0) : x(x), y(y), indices() {}

    Point(const Point &other) { *this = other; }

    ~Point() = default;

    Point &operator=(const Point &other) {
        if (this != &other) {
            this->x = other.x;
            this->y = other.y;
            this->indices = other.indices;
        }
        return *this;
    }

    string ToString() const {
        ostringstream convert;
        convert << "{";
        convert << "(" << this->x << "," << this->y << "),[";
        for (int i = 0; i < this->indices.size(); i++) {
            convert << this->indices[i];
            if (i != this->indices.size() - 1) {
                convert << ",";
            }
        }
        convert << "]";
        convert << "}";
        return convert.str();
    }

    friend bool operator==(const Point &lhs, const Point &rhs) {
        return (lhs.x == rhs.x) && (lhs.y == rhs.y);
    }

    friend bool operator!=(const Point &lhs, const Point &rhs) { return !(lhs == rhs); }

    friend Point operator+(const Point &p1, const Point &p2) { return Point(p1.x + p2.x, p1.y + p2.y); }

    friend Point operator-(const Point &p1, const Point &p2) { return Point(p1.x - p2.x, p1.y - p2.y); }

    friend ostream &operator<<(ostream &os, const Point &point) {
        os << point.ToString();
        return os;
    }
};

class Node {

public:
    float half_width = 0;
    float half_height = 0;
    Point center;
    Point bottom_left;
    Point top_right;

    explicit Node(const float x = 0, const float y = 0, const float half_width = 0, const float half_height = 0)
            : center(Point(x, y)), half_width(half_width), half_height(half_height) {
        // calculate bound
        Point half(half_width, half_height);
        this->bottom_left = center - half;
        this->top_right = center + half;
    }

    Node(const Node &other) { *this = other; }

    ~Node() = default;

    Node &operator=(const Node &other) {
        if (this != &other) {
            this->half_width = other.half_width;
            this->half_height = other.half_height;
            this->center = other.center;
            this->bottom_left = other.bottom_left;
            this->top_right = other.top_right;
        }
        return *this;
    }


    bool Contain(const float x, const float y) const {
        if (x <= bottom_left.x || x > top_right.x) {
            return false;
        }
        return !(y <= bottom_left.y || y > top_right.y);
    }

    bool Intersect(const Node &other) const {
        if (bottom_left.x > other.top_right.x) return false;
        if (other.bottom_left.x > top_right.x) return false;
        if (bottom_left.y > other.top_right.y) return false;
        return other.bottom_left.y <= top_right.y;
    }

    string ToString() const {
        ostringstream convert;
        convert << "Node(" << "Center - " << this->center << ", Width - " << this->half_width * 2 << ", Height - "
                << this->half_height * 2 << ")";
        return convert.str();
    }

    friend ostream &operator<<(ostream &os, const Node &quad) {
        os << quad.ToString();
        return os;
    }

    friend bool operator==(const Node &lhs, const Node &rhs) {
        return lhs.center == rhs.center && lhs.bottom_left == rhs.bottom_left && lhs.top_right == rhs.top_right;
    }

    friend bool operator!=(const Node &lhs, const Node &rhs) { return !(lhs == rhs); }

};


class QuadTree {
public:
    unique_ptr<Node> bounds;
    bool is_bound = false;
    vector<Point> nodes;
    unique_ptr<QuadTree> top_left;
    unique_ptr<QuadTree> top_right;
    unique_ptr<QuadTree> bottom_left;
    unique_ptr<QuadTree> bottom_right;

public:
    explicit QuadTree(const float x = 0, const float y = 0, const float half_width = 0, const float half_height = 0) {
        this->SetBound(x, y, half_width, half_height);
        this->Clear();
    }

    QuadTree(const QuadTree &other) { *this = other; }

    QuadTree(QuadTree &other) { *this = other; }

    QuadTree &operator=(const QuadTree &other) {
        if (this != &other) {
            this->bounds = make_unique<Node>(*other.bounds);
            this->nodes = other.nodes;
            if (other.top_left != nullptr) this->top_left = make_unique<QuadTree>(*other.top_left);
            if (other.top_right != nullptr) this->top_right = make_unique<QuadTree>(*other.top_right);
            if (other.bottom_left != nullptr) this->bottom_left = make_unique<QuadTree>(*other.bottom_left);
            if (other.bottom_right != nullptr) this->bottom_right = make_unique<QuadTree>(*other.bottom_right);
        }
        return *this;
    }


    bool Usable() const { return this->is_bound; }

    void SetBound(const float x, const float y, const float half_width, const float half_height) {
        this->bounds = make_unique<Node>(x, y, half_width, half_height);
        this->is_bound = true;
    }

    // insert a point
    bool Insert(const Point &p) {
        if (!this->bounds->Contain(p.x, p.y)) {
            return false;
        }

        if (this->nodes.size() < 4 && this->top_left == nullptr) {
            this->nodes.push_back(p);
            return true;
        }

        if (this->top_left == nullptr) {
            // divide the node
            const Point &center = this->bounds->center;
            float new_half_width = this->bounds->half_width / 2;
            float new_half_height = this->bounds->half_height / 2;
            float left = center.x - new_half_width;
            float right = center.x + new_half_width;
            float up = center.y + new_half_height;
            float down = center.y - new_half_height;

            this->top_left = make_unique<QuadTree>(left, up, new_half_width, new_half_height);
            this->top_right = make_unique<QuadTree>(right, up, new_half_width, new_half_height);
            this->bottom_left = make_unique<QuadTree>(left, down, new_half_width, new_half_height);
            this->bottom_right = make_unique<QuadTree>(right, down, new_half_width, new_half_height);

            for (Point &p : this->nodes) {
                if (this->top_left->Insert(p)) continue;
                if (this->top_right->Insert(p)) continue;
                if (this->bottom_left->Insert(p)) continue;
                if (this->bottom_right->Insert(p)) continue;
            }

            this->nodes.clear();
        }

        if (this->top_left->Insert(p)) return true;
        if (this->top_right->Insert(p)) return true;
        if (this->bottom_left->Insert(p)) return true;
        return this->bottom_right->Insert(p);
    }

    static float abs(float a) {
        return a >= 0 ? a : -a;
    }

    // recursively Query Point
    void Query(const float x, const float y, vector<const Point *> &output) {
        if (!this->bounds->Contain(x, y)) {
            return;
        }
        Point target(x, y);
        for (const Point &p : this->nodes) {
            Point diff = p - target;
            if (abs(diff.x) < 1e-4 && abs(diff.y) < 1e-4) {
                output.push_back(&p);
            }
        }
        if (this->top_left == nullptr) {
            return;
        }
        this->top_left->Query(x, y, output);
        this->top_right->Query(x, y, output);
        this->bottom_left->Query(x, y, output);
        this->bottom_right->Query(x, y, output);
    }

    // recursively Query Range
    void Query(const Node &range, vector<const Point *> &output) {
        if (!this->bounds->Intersect(range)) {
            return;
        }
        for (const Point &p : this->nodes) {
            if (range.Contain(p.x, p.y)) {
                output.push_back(&p);
            }
        }
        if (this->top_left == nullptr) {
            return;
        }
        this->top_left->Query(range, output);
        this->top_right->Query(range, output);
        this->bottom_left->Query(range, output);
        this->bottom_right->Query(range, output);
    }

    void Clear() {
        this->nodes.clear();
        this->top_left.reset();
        this->top_right.reset();
        this->bottom_left.reset();
        this->bottom_right.reset();
    }

    string ToString(const int level = 0) const {
        ostringstream convert;
        for (int i = 0; i < level; i++) {
            convert << "\t";
        }
        string tab = convert.str();
        convert.str("");
        convert.clear();
        convert << tab << "Bounds: " << *this->bounds << endl;
        convert << tab << "Points: [";
        for (int i = 0; i < this->nodes.size(); i++) {
            convert << this->nodes[i];
            if (i != this->nodes.size() - 1) {
                convert << ", ";
            }
        }
        convert << "]" << endl << endl;
        if (this->top_left != nullptr) {
            convert << this->top_left->ToString(level + 1);
        }
        if (this->top_left != nullptr) {
            convert << this->top_right->ToString(level + 1);
        }
        if (this->top_left != nullptr) {
            convert << this->bottom_left->ToString(level + 1);
        }
        if (this->top_left != nullptr) {
            convert << this->bottom_right->ToString(level + 1);
        }
        return convert.str();
    }

    friend ostream &operator<<(ostream &os, const QuadTree &tree) {
        os << tree.ToString();
        return os;
    }

};


template<typename K, typename V>
class HashMap {
public:
    const float rehash_threshold = 0.7;
    size_t longest_probe = 0;
    size_t capacity = 0;
    size_t size = 0;

    vector<pair<K, vector<V>>> hash_array;
    enum Status {
        EMPTY, OCCUPIED
    };
    vector<Status> status;


    size_t QuadraticProbing(const size_t hash, const size_t i) const {
        size_t n = hash + i;
        return ((n * n + n) / 2) % this->capacity;
    }


    void ReHash() {
        this->size = 0;
        this->capacity *= 2;
        vector<pair<K, vector<V>>> oldBuckets{move(this->hash_array)};
        vector<Status> oldStatus{move(this->status)};

        this->hash_array.resize(this->capacity);
        this->status.resize(this->capacity);

        for (int i = 0; i < oldBuckets.size(); i++) {
            if (oldStatus[i] == OCCUPIED) {
                const vector<V> &values = oldBuckets[i].second;
                for (const V value : values) {
                    this->Insert(oldBuckets[i].first, value);
                }
            }
        }
    }

    explicit HashMap() {

        this->capacity = 1024;
        this->hash_array.resize(this->capacity);
        this->status.resize(this->capacity);
    }

    HashMap(const HashMap &other) {
        *this = other;
    }

    ~HashMap() = default;

    HashMap &operator=(const HashMap &other) {
        if (&other != this) {
            this->capacity = other.capacity;
            this->size = other.size;
            this->hash_array.resize(this->capacity);
            this->status.resize(this->capacity);
            copy(other.hash_array.begin(), other.hash_array.end(), this->hash_array.begin());
            copy(other.status.begin(), other.status.end(), this->hash_array.begin());
        }
        return *this;
    }

    // Insert <key,value> to the hash map
    bool Insert(const K &key, const V &value) {
        // visited
        unordered_set<unsigned long> visited;

        if (this->size / this->capacity >= this->rehash_threshold) {
            this->ReHash();
        }

        size_t h = hash<K>{}(key) % this->capacity;
        size_t i = 0;
        size_t hi = h + i;

        while (status[hi] == OCCUPIED) {
            // the key already exists
            if (hash_array[hi].first == key) {
                hash_array[hi].second.push_back(value);
                return true;
            }
            visited.insert(hi);
            i++;
            // quadratic probing
            hi = this->QuadraticProbing(h, i);

            if (visited.find(hi) != visited.end()) {
                this->ReHash();
                visited.clear();
            }
        }
        // longest probe
        this->longest_probe = max(this->longest_probe, i);
        // hash array set
        hash_array[hi] = {key, vector<V>({value})};
        // set status
        status[hi] = OCCUPIED;
        this->size++;
        return true;
    }

    const vector<V> &Find(const K &key) const {
        unordered_set<unsigned long> visited;
        unsigned long h = hash<K>{}(key) % this->capacity;
        unsigned long i = 0;
        unsigned long hi = h + i;
        while (status[hi] != EMPTY) {
            if (status[hi] == OCCUPIED && hash_array[hi].first == key) {
                return hash_array[hi].second;
            }
            visited.insert(hi);
            i++;
            hi = this->QuadraticProbing(h, i);
            if (visited.find(hi) != visited.end()) {
                break;
            }
        }

        throw out_of_range("Key not exist.");
    }

    void Clear() {
        this->size = 0;
        this->capacity = 1024;
        this->hash_array.clear();
        this->hash_array.resize(this->capacity);
        this->status.clear();
        this->status.resize(this->capacity);
    }

    string ToString() const {
        ostringstream convert;
        convert << "Hash Map:" << endl;
        convert << "\tSize: " << this->size << ", Capacity: " << this->capacity << endl;
        for (int i = 0; i < status.size(); i++) {
            if (status[i] == OCCUPIED) {
                convert << "\t" << i << " : " << "{" << hash_array[i].first << ",[";
                const vector<V> &values = hash_array[i].second;
                for (int j = 0; j < values.size() - 1; j++) {
                    convert << values[j] << ",";
                }
                if (!values.empty()) {
                    convert << values.back();
                }
                convert << "]}" << endl;
            }
        }
        return convert.str();
    }

    friend ostream &operator<<(ostream &os, const HashMap<K, V> &map) {
        os << map.ToString();
        return os;
    }
};


class DMS {
public:
    string latitude;
    string longitude;

    DMS() {
        this->SetLatitude("");
        this->SetLongitude("");
    }

    DMS(const string &latitude, const string &longitude) {
        this->SetLatitude(latitude);
        this->SetLongitude(longitude);
    }

    void SetLatitude(const string &latitude) {
        if (latitude == "Unknown") {
            this->latitude = "";
        } else {
            this->latitude = latitude;
        }
    }

    void SetLongitude(const string &longitude) {
        if (longitude == "Unknown") {
            this->latitude = "";
        } else {
            this->longitude = longitude;
        }
    }

    bool Valid() const { return !this->latitude.empty() && !this->longitude.empty(); }

    string ToString() const {
        if (this->latitude.empty() && this->longitude.empty()) {
            return "|";
        }
        return this->latitude + "|" + this->longitude;
    }

    friend ostream &operator<<(ostream &os, const DMS &coord) {
        os << coord.ToString();
        return os;
    }
};

class DEC {
public:
    float latitude;
    float longitude;

public:
    DEC() {
        this->SetLatitude(0);
        this->SetLongitude(0);
    }


    DEC(const DMS &dms) {
        this->latitude = LatitudeDMSToDEC(dms.latitude);
        this->longitude = LongitudeDMSToDEC(dms.longitude);
    }

    void SetLatitude(const string &latitude) {
        if (latitude.empty()) {
            this->SetLatitude(0);
        } else {
            this->SetLatitude(stod(latitude));
        }
    }

    void SetLongitude(const string &longitude) {
        if (longitude.empty()) {
            this->SetLongitude(0);
        } else {
            this->SetLongitude(stod(longitude));
        }
    }

    void SetLatitude(const float latitude) { this->latitude = latitude; }

    void SetLongitude(const float longitude) { this->longitude = longitude; }


    bool Valid() const { return this->latitude != 0.0 && this->longitude != 0.0; }

    string ToString() const {
        if (this->latitude == 0.0 && this->longitude == 0.0) {
            return "|";
        }

        return to_string(latitude) + "|" + to_string(longitude);
    }

    friend ostream &operator<<(ostream &os, const DEC &coord) {
        os << coord.ToString();
        return os;
    }

    static float LatitudeDMSToDEC(const string &latitude) {
        float degrees;
        float minutes;
        float seconds;
        if (latitude.empty()) {
            return 0;
        } else {
            degrees = stod(latitude.substr(0, 2));
            minutes = stod(latitude.substr(2, 2)) / 60;
            seconds = stod(latitude.substr(4, 2)) / 3600;
            return (degrees + minutes + seconds) * (latitude.back() == 'N' ? 1 : -1);
        }
    }

    static float LongitudeDMSToDEC(const string &longitude) {
        float degrees;
        float minutes;
        float seconds;
        if (longitude.empty()) {
            return 0;
        } else {
            degrees = stod(longitude.substr(0, 3));
            minutes = stod(longitude.substr(3, 2)) / 60;
            seconds = stod(longitude.substr(5, 2)) / 3600;
            return (degrees + minutes + seconds) * (longitude.back() == 'E' ? 1 : -1);
        }
    }

    static float SecondsToDEC(const string coord) {
        return stod(coord) / 3600;
    }
};


class GeoDataRecord {
public:
    static const unordered_set<string> POP_TYPES;
    static const unordered_set<string> WATER_TYPES;
    static const unordered_set<string> STRUCTURE_TYPES;

    size_t offset;
    int feature_id = -1;
    string name_index;
    string name;
    string feature_class;
    string state_alpha;
    string state_numeric;
    string county_name;
    string county_numeric;
    DMS primary_dms;
    DEC primary_dec;
    DMS source_dms;
    DEC source_dec;
    unique_ptr<int> elevation_meters;
    unique_ptr<int> elevation_feet;
    string map_name;
    string date_created;
    string date_edited;

    explicit GeoDataRecord(const vector<string> &tokens, const size_t offset = -1) {
        if (tokens.size() < 19) {
            return;
        }

        this->offset = offset;
        this->feature_id = stoi(tokens[0]);
        this->name = tokens[1];
        this->feature_class = tokens[2];
        this->state_alpha = tokens[3];
        this->state_numeric = tokens[4];
        this->county_name = tokens[5];
        this->county_numeric = tokens[6];
        this->primary_dms.SetLatitude(tokens[7]);
        this->primary_dms.SetLongitude(tokens[8]);
        this->primary_dec.SetLatitude(tokens[9]);
        this->primary_dec.SetLongitude(tokens[10]);
        this->source_dms.SetLatitude(tokens[11]);
        this->source_dms.SetLongitude(tokens[12]);
        this->source_dec.SetLatitude(tokens[13]);
        this->source_dec.SetLongitude(tokens[14]);

        if (tokens[15].empty()) {
            this->elevation_meters.reset(nullptr);
        } else {
            this->elevation_meters = make_unique<int>(stoi(tokens[15]));
        }

        if (tokens[16].empty()) {
            this->elevation_feet.reset(nullptr);
        } else {
            this->elevation_feet = make_unique<int>(stoi(tokens[16]));
        }

        this->map_name = tokens[17];
        this->date_created = tokens[18];

        if (tokens.size() == 20) {
            this->date_edited = tokens[19];
        }

        // generate name index
        ostringstream convert;
        convert << this->name << "|" << this->state_alpha;
        this->name_index = convert.str();

        // create valid DEC coordinates
        if (!this->primary_dec.Valid()) {
            this->primary_dec = DEC(this->primary_dms);
        }

        if (!this->source_dec.Valid()) {
            this->source_dec = DEC(this->source_dms);
        }
    }

    GeoDataRecord(const GeoDataRecord &other) { *this = other; }


    GeoDataRecord &operator=(const GeoDataRecord &other) {
        if (this != &other) {
            this->offset = other.offset;
            this->feature_id = other.feature_id;
            this->name_index = other.name_index;
            this->name = other.name;
            this->feature_class = other.feature_class;
            this->state_alpha = other.state_alpha;
            this->state_numeric = other.state_numeric;
            this->county_name = other.county_name;
            this->county_numeric = other.county_numeric;
            this->primary_dms = other.primary_dms;
            this->primary_dec = other.primary_dec;
            this->source_dms = other.source_dms;
            this->source_dec = other.source_dec;
            if (other.elevation_meters) {
                this->elevation_meters = make_unique<int>(*other.elevation_meters);
            }
            if (other.elevation_feet) {
                this->elevation_feet = make_unique<int>(*other.elevation_feet);
            }
            this->map_name = other.map_name;
            this->date_created = other.date_created;
            this->date_edited = other.date_edited;
        }
        return *this;
    }


    string ToString() const {
        ostringstream convert;
        convert << to_string(this->feature_id) << "|" << this->name << "|" << this->feature_class << "|";
        convert << this->state_alpha << "|" << this->state_numeric << "|";
        convert << this->county_name << "|" << this->county_numeric << "|";
        convert << this->primary_dms << "|" << this->primary_dec << "|";
        convert << this->source_dms << "|" << this->source_dec << "|";
        if (this->elevation_meters != nullptr) { convert << to_string(*this->elevation_meters); }
        convert << "|";
        if (this->elevation_feet != nullptr) { convert << to_string(*this->elevation_feet); }
        convert << "|";
        convert << this->map_name << "|";
        convert << this->date_created << "|";
        convert << this->date_edited << "|";
        return convert.str();
    }

    string ToLongString() const {
        ostringstream convert;

        if (this->feature_id != -1) {
            convert << "Feature ID: " << this->feature_id << endl;
        }

        if (!this->name.empty()) {
            convert << "Feature Name: " << this->name << endl;
        }

        if (!this->feature_class.empty()) {
            convert << "Feature Class: " << this->feature_class << endl;
        }

        if (!this->state_alpha.empty()) {
            convert << "State Alpha: " << this->state_alpha << endl;
        }

        if (!this->state_numeric.empty()) {
            convert << "State Numeric: " << this->state_numeric << endl;
        }

        if (!this->county_name.empty()) {
            convert << "County Name: " << this->county_name << endl;
        }

        if (!this->county_numeric.empty()) {
            convert << "County Numeric: " << this->county_numeric << endl;
        }

        if (this->primary_dms.Valid()) {
            convert << "Primary Coordinate(DMS): " << this->primary_dms << endl;
        }

        if (this->primary_dec.Valid()) {
            convert << "Primary Coordinate(DEC): " << this->primary_dec << endl;
        }

        if (this->source_dms.Valid()) {
            convert << "Source Coordinate(DMS): " << this->source_dms << endl;
        }

        if (this->source_dec.Valid()) {
            convert << "Source Coordinate(DEC): " << this->source_dec << endl;
        }

        if (this->elevation_meters != nullptr) {
            convert << "Elevation in Meters: " << *this->elevation_meters << endl;
        }

        if (this->elevation_feet != nullptr) {
            convert << "Elevation in Feet: " << *this->elevation_feet << endl;
        }

        if (!this->map_name.empty()) {
            convert << "Map Name: " << this->map_name << endl;
        }

        if (!this->date_created.empty()) {
            convert << "Date Created: " << this->date_created << endl;
        }

        if (!this->date_edited.empty()) {
            convert << "Date Edited: " << this->date_edited << endl;
        }

        return convert.str();
    }

    friend ostream &operator<<(ostream &os, const GeoDataRecord &gf) {
        os << gf.ToString();
        return os;
    }
};


const unordered_set<string> GeoDataRecord::POP_TYPES = {
        "Populated Place"
};

const unordered_set<string> GeoDataRecord::WATER_TYPES = {
        "Arroyo", "Basin", "Bay", "Bend", "Canal", "Channel", "Falls",
        "Glacier", "Gut", "Harbor", "Lake", "Rapids", "Reservoir",
        "Sea", "Spring", "Stream", "Swamp", "Well"
};

const unordered_set<string> GeoDataRecord::STRUCTURE_TYPES = {
        "Airport", "Bridge", "Building", "Church", "Dam", "Hospital",
        "Levee", "Park", "Post Office", "School", "Tower", "Tunnel"
};


// Cache Item
class CacheItem {
public:
    chrono::nanoseconds timestamp;
    const GeoDataRecord feature;

    CacheItem(const GeoDataRecord &gdf) : feature(gdf) {
        this->UpdateTimestamp();
    }

    void UpdateTimestamp() {
        this->timestamp = chrono::duration_cast<chrono::nanoseconds>(
                chrono::system_clock::now().time_since_epoch());
    }

    string ToString() const {
        ostringstream convert;
        convert << "{";
        convert << "timestamp:" << chrono::duration_cast<chrono::seconds>(this->timestamp).count() << ",";
        convert << "offset:" << this->feature.offset << ",";
        convert << "feature_id:" << this->feature.feature_id << ",";
        convert << "name_index:" << this->feature.name_index;
        convert << "}";
        return convert.str();
    }

    friend ostream &operator<<(ostream &os, const CacheItem &ce) {
        os << ce.ToString();
        return os;
    }
};


class GeoDataBase {
public:

    static constexpr int CACHE_SIZE = 15;
    list<CacheItem> cache;
    fstream out_file;
    HashMap<string, size_t> name_index;
    QuadTree coordinate_index;

    explicit GeoDataBase(const string &fn) {
        if (fn.empty()) {
            return;
        }
        // clear already opened
        if (this->out_file.is_open()) {
            this->out_file.close();
            this->cache.clear();
            this->name_index.Clear();
            this->coordinate_index.Clear();
        }

        this->out_file.open(fn, fstream::in | fstream::out | fstream::trunc);
    }

    ~GeoDataBase() {
        if (this->out_file.is_open()) {
            this->out_file.close();
        }
    }

    size_t ReadData(const string &filename) {
        size_t count = 0;
        string line;
        ifstream file;
        file.open(filename);
        if (!file.is_open()) {
            return count;
        }
        // skip header
        getline(file, line);
        while (!file.eof()) {
            getline(file, line);
            if (line.empty()) { continue; }
            this->WriteLine(line);
            count++;
        }
        file.close();
        return count;
    }


    bool WriteLine(const string &line) {
        auto entry = GeoDataRecord(Split(line, "|"));
        if (!this->coordinate_index.Usable()) {
            return false;
        }
        size_t init_pos = this->out_file.tellg();
        this->out_file << entry << endl;
        // insert hash map
        name_index.Insert(entry.name_index, init_pos);
        // insert quad tree
        Point p(entry.primary_dec.longitude, entry.primary_dec.latitude);
        p.indices.push_back(init_pos);
        this->coordinate_index.Insert(p);
        return true;
    }


    vector<GeoDataRecord> NameSearch(const string &name, const string &state) {
        return this->Retrive(this->name_index.Find(name + "|" + state));
    }

    vector<GeoDataRecord> CoordinateSearch(const DMS &coord) {
        return this->CoordinateSearch(DEC(coord));
    }

    vector<GeoDataRecord> CoordinateSearch(const DEC &coord) {
        vector<GeoDataRecord> features;
        if (!this->coordinate_index.Usable()) {
            return features;
        }
        vector<const Point *> points;
        this->coordinate_index.Query(coord.longitude, coord.latitude, points);
        this->PointsToGeoFeatures(points, features);
        return features;
    }

    vector<GeoDataRecord> CoordinateSearch(const DMS &coord, const float half_width, const float half_height) {
        return this->CoordinateSearch(DEC(coord), half_width, half_height);
    }

    vector<GeoDataRecord> CoordinateSearch(const DEC &coord, const float half_width, const float half_height) {
        vector<GeoDataRecord> features;
        if (!this->coordinate_index.Usable()) {
            return features;
        }
        Node range(coord.longitude, coord.latitude, half_width, half_height);
        vector<const Point *> points;
        this->coordinate_index.Query(range, points);
        this->PointsToGeoFeatures(points, features);
        return features;
    }

    string ToString() const {
        int i = 1;
        ostringstream convert;
        convert << "Buffer Pool:" << endl;
        convert << "MRU" << endl;
        for (const CacheItem &ce : this->cache) {
            convert << i++ << ": " << ce << endl;
        }
        convert << "LRU" << endl;
        return convert.str();
    }


    vector<GeoDataRecord> Retrive(const vector<size_t> &offsets) {
        string line;
        vector<GeoDataRecord> features;

        for (const size_t offset : offsets) {
            // test if data in cache
            bool in_cache = false;
            for (auto it = this->cache.begin(); it != this->cache.end(); ++it) {
                if (it->feature.offset == offset) {
                    features.push_back(it->feature);
                    // update its last use time
                    it->UpdateTimestamp();
                    this->cache.sort([](
                            const CacheItem &a, const CacheItem &b) {
                        return a.timestamp > b.timestamp;
                    });
                    in_cache = true;
                }
            }
            // already in cache
            if (in_cache) {
                continue;
            }

            // not in cache , retrive from file
            this->out_file.seekg(offset);
            getline(this->out_file, line);

            const GeoDataRecord result(Split(line, "|"), offset);
            this->UpdateCache(result);

            features.push_back(result);
        }

        return features;
    }

    void UpdateCache(const GeoDataRecord &entry) {
        if (this->cache.size() >= CACHE_SIZE) {
            this->cache.pop_back();
        }

        this->cache.emplace_front(entry);
        this->cache.sort([](
                const CacheItem &a, const CacheItem &b) {
            return a.timestamp > b.timestamp;
        });
    }

    void PointsToGeoFeatures(const vector<const Point *> &points, vector<GeoDataRecord> &output) {
        set<size_t> offsets;
        for (const Point *p : points) {
            for (const size_t i : p->indices) {
                offsets.insert(i);
            }
        }
        vector<GeoDataRecord> features = this->Retrive(
                vector<size_t>(offsets.begin(), offsets.end()));

        output.clear();
        output.insert(output.end(), features.begin(), features.end());
    }
};

// Command handler
class Commands {
public:
    const unordered_map<string, int> NUM_OF_ARGS = {
            {"world",      4},
            {"import",     1},
            {"debug",      1},
            {"quit",       0},
            {"what_is_at", 2},
            {"what_is",    2},
            {"what_is_in", 4},
            {";",          0}
    };
    string cmd = "";
    vector<string> args;

public:
    explicit Commands(const vector<string> &tokens) {
        if (tokens.empty()) {
            return;
        }
        if (tokens[0][0] == ';') {
            cmd = ";";
            args.push_back(tokens[0]);
            return;
        }
        if (NUM_OF_ARGS.find(tokens[0]) == NUM_OF_ARGS.end()) {
            return;
        }
        cmd = tokens[0];
        if (tokens.size() == 1) {
            return;
        }
        for (int i = 1; i < tokens.size(); i++) {
            this->args.push_back(tokens[i]);
        }
    }

    bool Valid() const {
        if (NUM_OF_ARGS.find(this->cmd) == NUM_OF_ARGS.end()) { return false; }
        if (this->cmd == ";") {
            return true;
        } else {
            return this->args.size() >= NUM_OF_ARGS.at(this->cmd);
        }
    }

    string ToString() const {
        ostringstream convert;
        if (this->cmd == ";") {
            for (const string &s : this->args) {
                convert << s;
            }
        } else {
            convert << this->cmd << "\t";
            for (const auto &arg : this->args) {
                convert << arg << "\t";
            }
        }
        return convert.str();
    }

    friend ostream &operator<<(ostream &os, const Commands &cmd) {
        os << cmd.ToString();
        return os;
    }
};

// Line Tokenizer
class Tokenizer {
public:
    string delimiter;
    string filename;
    ifstream input_file;

public:
    Tokenizer(const string &filename, const string &delimiter)
            : delimiter(delimiter), filename(filename), input_file(filename) {}

    ~Tokenizer() {
        this->input_file.close();
    }

    vector<string> LineToken(const bool skip_empty = false) {
        if (!this->input_file.is_open()) {
            return vector<string>();
        }
        string line;
        getline(this->input_file, line);
        if (line[0] == ';') {
            return {line};
        }
        return Split(line, this->delimiter, skip_empty);
    }
};


class GIS {
public:
    vector<Commands> command_vector;
    GeoDataBase db;
    fstream log_file;
    int command_state = -1;

    GIS(const string &database_file, const string &script_file, const string &log_file) : db(database_file) {
        // open log file
        if (this->log_file.is_open()) {
            this->log_file.close();
        }
        this->log_file.open(log_file, fstream::out | fstream::trunc);
        this->ParseCommand(script_file);
        
        ostringstream convert;
        convert << "GeoDataBase File: " << database_file << endl;
        convert << "Command Script: " << script_file << endl;
        convert << "Log File: " << log_file << endl;
        this->log_file << convert.str() << endl;
        // execute command
        for (const Commands &cmd : this->command_vector) {
            this->Execute(cmd);
        }
    }

    ~GIS() = default;
    
    void ParseCommand(const string &script) {
        Tokenizer command_tokenizer(script, "\t");
        while (!command_tokenizer.input_file.eof()) {
            Commands cmd(command_tokenizer.LineToken());
            if (cmd.Valid()) {
                this->command_vector.push_back(cmd);
            }
        }
    }
    
    void Execute(const Commands &cmd) {
        if (!cmd.Valid()) { return ; }

        const string &instruction = cmd.cmd;
        const vector<string> &args = cmd.args;

        this->LogCommand(cmd);

        if (instruction == ";") {
            return ;
        }

        if (instruction == "quit") {
            return ;
        }
        auto startTime = chrono::high_resolution_clock::now();
        // command
        if (instruction == "world") {
            this->NewWorld(args[0], args[1], args[2], args[3]);
        } else if (instruction == "import") {
            this->ImportRecord(args[0]);
        } else if (instruction == "debug") {
            this->Debug(args[0]);
        } else if (instruction == "what_is_at") {
            this->CoordinateSearch(args[0], args[1]);
        } else if (instruction == "what_is") {
            this->SearchName(args[0], args[1]);
        } else if (instruction == "what_is_in") {
            if (args[0] == "-long") {
                this->SearchNode(args[1], args[2], args[3], args[4], true);
            } else if (args[0] == "-filter") {
                this->SearchNode(args[2], args[3], args[4], args[5], false, args[1]);
            } else {
                this->SearchNode(args[0], args[1], args[2], args[3]);
            }
        }
        auto endTime = chrono::high_resolution_clock::now();
        chrono::duration<double> delta = endTime - startTime;
        double time = delta.count();
        log_file << ("Time elapsed: " + to_string(time) + "s") << endl;
        log_file << ("-----------------------------------------------------------------") << endl;
        
    }

    bool NewWorld(const string &west, const string &east, const string &south, const string &north) {
        const float dec_west = DEC::LongitudeDMSToDEC(west);
        const float dec_east = DEC::LongitudeDMSToDEC(east);
        const float dec_south = DEC::LatitudeDMSToDEC(south);
        const float dec_north = DEC::LatitudeDMSToDEC(north);

        const float half_width = dec_east - dec_west;
        const float half_height = dec_north - dec_south;

        const float center_x = dec_west + half_width;
        const float center_y = dec_south + half_height;

        this->db.coordinate_index.SetBound(center_x, center_y, half_width, half_height);
        ostringstream convert;
    convert << "World boundaries are set to:"<<endl;
    convert << "           " << north <<"              " << endl;
    convert << west << "                    " << east << endl;
    convert <<"           " << south <<"              " << endl;
        log_file << convert.str() << endl;
        return true;
    }


    bool ImportRecord(const string &filename) {
        size_t results = this->db.ReadData(filename);
        ostringstream convert;
        convert << "Imported " << results << " records from " << filename << endl;
        convert << "Longest Probe: " << this->db.name_index.longest_probe << endl;
        convert << "Imported Locations: " << results << endl;
        log_file << convert.str() << endl;
        return true;
    }

    bool Debug(const string &option) {
        if (option == "quad") {
            this->DebugQuadTree();
        } else if (option == "hash") {
            this->DebugHashMap();
        } else if (option == "pool") {
            this->DebugPool();
        } else {
            return false;
        }
        return true;
    }

    bool CoordinateSearch(const string &latitude, const string &longitude) {
        vector<GeoDataRecord> features = this->db.CoordinateSearch(DMS(latitude, longitude));
        sort(features.begin(), features.end(), [](auto &a, auto &b) { return a.name < b.name; });
        ostringstream convert;
        convert << "Results: " << features.size() << endl;
        for (const GeoDataRecord &feature : features) {
            convert << feature.offset << " " << feature.name << " " << feature.county_name << " "
                    << feature.state_alpha << endl;
        }
        log_file << convert.str() << endl;
        return true;
    }


    bool SearchName(const string &name, const string &state) {
        try {
            vector<GeoDataRecord> features = this->db.NameSearch(name, state);
            sort(features.begin(), features.end(), [](auto &a, auto &b) { return a.name < b.name; });
            ostringstream convert;
            convert << "Results: " << features.size() << endl;
            for (const GeoDataRecord &feature : features) {
                convert << to_string(feature.offset) << " " << feature.county_name << " "
                        << feature.primary_dms << endl;
            }
            log_file << convert.str() << endl;
            return true;
        } catch (const exception &e) {
            return false;
        }
    }

    bool SearchNode(const string &latitude, const string &longitude,
                      const string &halfLat, const string &halfLng,
                      const bool longFormat = false, const string &filter = "") {

        vector<GeoDataRecord> output;
        const vector<GeoDataRecord> &features =
                this->db.CoordinateSearch(DMS(latitude, longitude), DEC::SecondsToDEC(halfLng),
                                          DEC::SecondsToDEC(halfLat));
        // populated type
        if (filter == "pop") {
            copy_if(features.begin(), features.end(), back_inserter(output), [](const GeoDataRecord &f) {
                return GeoDataRecord::POP_TYPES.find(f.feature_class) != GeoDataRecord::POP_TYPES.end();
            });
        }
        // water type
        else if (filter == "water") {
            copy_if(features.begin(), features.end(), back_inserter(output), [](const GeoDataRecord &f) {
                return GeoDataRecord::WATER_TYPES.find(f.feature_class) != GeoDataRecord::WATER_TYPES.end();
            });
        }
        // structure type
        else if (filter == "structure") {
            copy_if(features.begin(), features.end(), back_inserter(output), [](const GeoDataRecord &f) {
                return GeoDataRecord::STRUCTURE_TYPES.find(f.feature_class) != GeoDataRecord::STRUCTURE_TYPES.end();
            });
        }
        // other type
        else {
            output = features;
        }


        sort(output.begin(), output.end(), [](auto &a, auto &b) { return a.name < b.name; });

        ostringstream convert;

        convert << "Results: " << output.size() << endl;

        for (const GeoDataRecord &feature : output) {
            if (longFormat) {
                convert << feature.ToLongString() << endl;
            } else {
                convert << to_string(feature.offset) << " "
                        << feature.name << " "
                        << feature.state_alpha << " "
                        << feature.primary_dms << endl;
            }
        }
        log_file << convert.str() << endl;

        return true;
    }


    void LogCommand(const Commands &command) {
        ostringstream convert;
        if (command.cmd == ";") {
            convert << command;
        } else {
            this->command_state++;
            convert << "Command " << this->command_state << ": " << command << endl;
        }
        log_file << convert.str() << endl;
    }
    void DebugQuadTree() {
        log_file << this->db.coordinate_index.ToString() << endl;
    }
    void DebugHashMap() {
        log_file << this->db.name_index.ToString() << endl;
    }
    void DebugPool() {
        log_file << this->db.ToString() << endl;
    }
};


int main(int argc, char *argv[]) {
    if (argc != 4) {
        cout << "Usage: ./" << argv[0] << " <database file name> <command script file name> <log file name>" << endl;
        exit(-1);
    }
    GIS gis(argv[1], argv[2], argv[3]);
    return 0;
}
