#include <chrono>

namespace timer
{
    namespace time = std::chrono;

    class timer_t
    {
      public:
        timer_t() : m_start(time::system_clock::now()) {}

        template<typename T>
        T get_elapsed_sec()
        {
            return time::duration_cast<time::duration<T>>(
                       time::system_clock::now() - m_start)
                .count();
        }

      private:
        time::time_point<time::system_clock> m_start;
    };
}    // namespace timer
