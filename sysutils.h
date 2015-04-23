#ifndef SYSUTILS
#define SYSUTILS

#include <string>
#include <vector>
#include <dirent.h>
#include <time.h>
#include <sys/stat.h>


namespace core
{
    template <class T>
    std::string to_string(T value)
    {
        std::stringstream ss;
        ss << value;
        std::string str;
        ss >> str;
        return str;
    }

    struct Time
    {
        short sec;
        short min;
        short hour;
        short day;
        short month;
        short year;

        Time()  {}

        std::string timeString()
        {
            return to_string(hour) + "." +
                   to_string(min) + "." +
                   to_string(sec);
        }

        void operator++()
        {
            if (++sec > 59)
            {
                sec = 0;
                if (++min > 59)
                {
                    min = 0;
                    if (++hour > 23)
                        hour = 0;
                }

            }
        }
    };

    std::string getExt(const std::string& path)
    {
        size_t pose = path.find_last_of(".");
        if (std::string::npos == pose)
            return path;
        return std::string(path.begin() + pose + 1, path.end());
    }

    bool isFileExists(const char* path)
    {
        struct stat buffer;
        return (stat (path, &buffer) == 0); 
    }

    bool isDir(const char* path)
    {
        struct stat sb;
        return stat(path, &sb) == 0 && S_ISDIR(sb.st_mode);
    }

    bool mkDir(const char* path)
    {
        return !isDir(path) ?
            mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) :
            false;
    }

    bool listFiles(const std::string& path, std::vector<std::string>& files)
    {
        files.clear();
        DIR *dir;
        struct dirent *ent;
        if ((dir = opendir (path.c_str())) != NULL)
        {
            /* print all the files and directories within directory */
            while ((ent = readdir (dir)) != NULL)
            {
                if (strncmp(ent->d_name, ".", 1) == 0 ||
                    strncmp(ent->d_name, "..", 2) == 0)
                    continue;
                files.push_back(std::string(ent->d_name));
            }
            closedir (dir);
            return true;
        }
        else
        {
          /* could not open directory */
          perror ("");
          return false;
        }
    }

    std::string nameWithExt(const std::string& path)
    {
        size_t pose = path.find_last_of("\\/");
        if (std::string::npos == pose)
            return path;
        return std::string(path.begin() + pose + 1, path.end());
    }

    std::string name(const std::string& path)
    {
        size_t pose1 = path.find_last_of("\\/");
        size_t pose2 = path.find_last_of(".");
        if (std::string::npos == pose1)
            return (std::string::npos == pose2) ? path : std::string(path.begin(), path.begin() + pose2);

        return std::string(path.begin() + pose1 + 1, path.begin() + pose2);
    }

    void getTime(Time& t)
    {
        time_t rawtime;
        struct tm  *timeinfo;
    
        time (&rawtime);
        timeinfo = localtime(&rawtime);
    
        t.sec = timeinfo->tm_sec;
        t.min = timeinfo->tm_min;
        t.hour = timeinfo->tm_hour;
    
        t.day = timeinfo->tm_mday;
        t.month = timeinfo->tm_mon + 1;
        t.year = 1900 + timeinfo->tm_year;
    }
    
    std::string getTime()
    {
        time_t rawtime;
        struct tm  *timeinfo;
        char buffer[80];
    
        time (&rawtime);
        timeinfo = localtime(&rawtime);
    
        //strftime(buffer, 80, "%I.%M.%S", timeinfo);
        strftime(buffer, 80, "%d-%m-%Y %I:%M:%S", timeinfo);
        //std::string str = std::string(timeinfo->tm_hour) + "." + std::string(timeinfo->tm_min) + "." + std::string(timeinfo->tm_sec);
        std::string str(buffer);
    
        std::string newstr(str.begin() + 11, str.end());
        int pos = newstr.find_last_of(":");
        newstr[pos] = '.';
        pos = newstr.find_last_of(":");
        newstr[pos] = '.';
        //std::cout << newstr << std::endl;
        return newstr;
    
        //return str;
    }
};

#endif // SYSUTILS

