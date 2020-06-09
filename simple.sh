ualpha:
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 10 -o .\result\ualpha.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 409600 -s ualpha

smis_jacobian
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o .\result\smis_jacobian_k2_stratified.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s smis_jacobian -k 2 -x

smis_all
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o .\result\smis_all_k2_stratified.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s smis_all -k 4 -x

cmis
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o .\result\cmis.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s cmis

average
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o .\result\Average.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s average

proxy
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o .\result\proxy.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s proxy_sample

proxy-sum
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o .\result\proxy_sum.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s proxy_sample -x

cargo run --features="pbrt progress-bar" --release --example=cli -- -t -1 -n 1 -o .\result\proxy.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 1 -s proxy_sample >a.txt

sims Elapsed Integrator: 27328 ms 2
smis Elapsed Integrator: 31123 ms 3


cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o .\result\smis_all_k2_stratified.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s smis_all -k 2 -x


cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o .\result\proxy_1.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 10240 -s proxy_sample;
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 2 -o .\result\proxy_2.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 10240 -s proxy_sample;
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 3 -o .\result\proxy_3.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 10240 -s proxy_sample;
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 4 -o .\result\proxy_4.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 10240 -s proxy_sample;
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 5 -o .\result\proxy_5.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 10240 -s proxy_sample;
$echo