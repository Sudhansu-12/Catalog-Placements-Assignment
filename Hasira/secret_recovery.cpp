#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include "json.hpp"
#include <gmpxx.h>

using json = nlohmann::json;

mpz_class decode_base(const std::string &num_str, int base)
{
    static const std::string digits = "0123456789abcdefghijklmnopqrstuvwxyz";
    static std::map<char, int> digit_map;

    if (digit_map.empty())
    {
        for (size_t i = 0; i < digits.size(); ++i)
        {
            digit_map[digits[i]] = i;
        }
    }

    mpz_class result = 0;
    for (char c : num_str)
    {
        char clower = std::tolower(static_cast<unsigned char>(c));
        if (digit_map.find(clower) == digit_map.end())
        {
            throw std::runtime_error(std::string("Invalid digit '") + c + "' in number");
        }
        int val = digit_map[clower];
        if (val >= base)
        {
            throw std::runtime_error("Digit out of base range");
        }
        result = result * base + val;
    }
    return result;
}

std::vector<std::pair<mpz_class, mpz_class>> get_points(const json &testcase)
{
    int n = testcase.at("keys").at("n").get<int>();
    int k = testcase.at("keys").at("k").get<int>();

    std::vector<std::pair<mpz_class, mpz_class>> points;
    for (auto it = testcase.begin(); it != testcase.end(); ++it)
    {
        if (it.key() == "keys")
            continue;

        const std::string &key = it.key();
        mpz_class x(key);

        const json &obj = it.value();
        int base = std::stoi(obj.at("base").get<std::string>());
        std::string val_str = obj.at("value");

        mpz_class y = decode_base(val_str, base);

        points.emplace_back(x, y);
    }

    std::sort(points.begin(), points.end(),
              [](const auto &a, const auto &b)
              { return a.first < b.first; });

    if ((int)points.size() < k)
        throw std::runtime_error("Not enough points for interpolation");

    points.resize(k);
    return points;
}

mpz_class lagrange_interpolate_constant(const std::vector<std::pair<mpz_class, mpz_class>> &points)
{
    int k = points.size();
    mpq_class secret = 0;

    for (int j = 0; j < k; ++j)
    {
        mpq_class numerator = 1;
        mpq_class denominator = 1;

        for (int m = 0; m < k; ++m)
        {
            if (m == j)
                continue;
            numerator *= mpq_class(-points[m].first);
            denominator *= mpq_class(points[j].first - points[m].first);
        }

        mpq_class lj = numerator / denominator;
        secret += mpq_class(points[j].second) * lj;
    }

    secret.canonicalize();

    if (secret.get_den() != 1)
    {
        std::cerr << "Warning: Result is fractional!\n";
    }

    return secret.get_num();
}

int main()
{
    try
    {
        std::ifstream infile("input.json");
        if (!infile.is_open())
        {
            std::cerr << "Could not open input.json\n";
            return 1;
        }

        json root;
        infile >> root;

        if (!root.is_array())
        {
            std::cerr << "JSON input must be an array of testcases\n";
            return 1;
        }

        for (const auto &testcase : root)
        {
            auto points = get_points(testcase);
            mpz_class secret = lagrange_interpolate_constant(points);
            std::cout << secret.get_str() << "\n";
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
