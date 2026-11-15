library(xfun)
litedown::reactor(attr.source = 'r')

z = strict_list(aaa = "I am aaa", b = 1:5)
z$a    # NULL — no partial match
z$aaa  # correct element

z2 = unclass(z)  # ordinary list
z2$a   # "I am aaa" — partial match fires!

z3 = as_strict_list(z2)
z3$a   # NULL again

x = structure(list(), foobar = 42)
base::attr(x, "foo")  # 42 — partial match
attr2(x, "foo")       # NULL
attr2(x, "foobar")    # 42

x = c("line one", "say \"hello\"", "tab\there")
x            # with indices and escapes
raw_string(x)  # what the strings actually contain

f = system.file("LICENSE", package = "xfun")
file_string(f)

fit = lm(dist ~ speed, data = cars)
tree(fit)

format_bytes(c(500, 2048, 1.5e6, 3e9))

grep_sub('(\\w+)@(\\w+)', '\\2/\\1', c('alice@corp', 'bob@org', 'not-an-email'))

is_blank(c("", "   ", "\t\n", "hello"))

join_words(c("apples", "oranges", "pears"))
join_words(c("R", "Python", "Julia"), and = "or")
join_words(c("one", "two", "three"), before = "'", after = "'")

alnum_id(c("Hello world 123!", "a  &b*^c"))
alnum_id("<h1>Section <em>1.2</em></h1>")

valid_syntax("1 + 1")
valid_syntax("1 +")       # incomplete
valid_syntax("if(T) {1 + 1")  # unclosed brace

code = c("# comment", "if (TRUE) {", "  1 + 1", "}", "print('done')")
split_source(code)
split_source(code, merge_comments = TRUE)

n2w(c(0, 1, 13, 42, 100, 1001), cap = TRUE)
n2w(1e6 + 234567, and = TRUE)
n2w(-987654321)

html_tag("a", "R Project", href = "https://www.r-project.org", target = "_blank")
html_tag("br")
html_tag("p", list("Visit ", html_tag("a", "CRAN", href = "https://cran.r-project.org")))

html_escape('2 < 3 & "quotes"')
html_escape('attr="val"', attr = TRUE)

strip_html('<p>Hello <strong>world</strong><!-- comment -->!</p>')

f = tempfile()
write_utf8(c("héllo", "wörld"), f)
read_utf8(f)

f = tempfile()
write_utf8("apple", f)
append_unique(c("apple", "banana", "cherry"), f)
read_utf8(f)

scripts = system.file("scripts", c("call-fun.R", "child-pids.sh"), package = "xfun")
read_all(scripts, before = function(f) paste("##", f), after = "")

f = tempfile()
write_utf8(c("color", "favourite", "grey"), f)
gsub_file(f, "colour|grey", c(colour = "color", grey = "gray"))
read_utf8(f)

f = tempfile()
write_utf8(c("a", "b", "c"), f)
process_file(f, rev)  # reverse the order of lines
read_utf8(f)

p = c("report.Rmd", "archive.tar.gz", "notes.nb.html", "script.R")
file_ext(p)
sans_ext(p)
with_ext(p, ".html")
with_ext(p, c(".pdf", ".zip", ".html", ".txt"))

normalize_path("~")
same_path("~/foo", file.path(Sys.getenv("HOME"), "foo"))

raw_bytes = charToRaw("Hello, xfun!")
(b64 = base64_encode(raw_bytes))
rawToChar(base64_decode(b64))

# data URI for a small file
f = system.file("LICENSE", package = "xfun")
substr(base64_uri(f), 1, 60)

tojson(list(x = 1:3, y = c(TRUE, FALSE), z = NULL))
tojson(head(iris, 3))    # named data frame → columns
tojson(matrix(1:6, 2))  # matrix → rows

tojson(list(n = 42, callback = js("function(x) { return x * 2; }")))

tojson(list(d = Sys.Date(), x = c(1, Inf, -Inf)), mode = "json")

x = rep(c("apple", "banana", "cherry"), 100)
nchar(tojson(x, pretty = FALSE))          # without dict encoding
nchar(tojson(x, dict = 1, pretty = FALSE))  # with dict encoding

path = tempfile()
x = 42
result = cache_exec(
  { Sys.sleep(0); x^2 },
  path = path, id = "square"
)
result  # 1764

# change x — cache is invalidated
x = 10
result = cache_exec(
  { Sys.sleep(0); x^2 },
  path = path, id = "square"
)
result  # 100

code = c(
  "1 + 1",
  "message('hello from record()')",
  "warning('a warning')",
  "head(mtcars, 2)"
)
res = record(code, dev.path = FALSE, error = TRUE)
tree(res)

pkg_attach("stats", "utils")  # equivalent to two library() calls

loadable("stats")
loadable("nonexistent.package.xyz")

Sys.setenv(R_XFUN_TEST_OPTION = "from-env")
env_option("xfun.test.option")       # "from-env"

options(xfun.test.option = TRUE)
env_option("xfun.test.option")       # TRUE (options() wins)
options(xfun.test.option = NULL)     # reset
env_option("xfun.test.option")       # "from-env" again
Sys.unsetenv("R_XFUN_TEST_OPTION")

old = set_envvar(c(MY_TOKEN = "secret"))
Sys.getenv("MY_TOKEN")
set_envvar(old)  # restore
Sys.getenv("MY_TOKEN")

in_dir(tempdir(), {
  writeLines("hello", "tmp.txt")
  file.exists("tmp.txt")
})
file.exists(file.path(tempdir(), "tmp.txt"))

is_windows()
is_macos()
is_linux()
is_unix()   # TRUE on macOS and Linux
is_arm64()  # TRUE on Apple silicon

try_error(log(1))
try_error(log("a"))

z = try_silent(sqrt("x"))
inherits(z, "try-error")

session_info(c("xfun", "litedown"), dependencies = FALSE)

for (i in 1:3) {
  do_once(message("This prints once."), "my-one-time-task")
}

