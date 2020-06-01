ualpha:
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o ualpha.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s ualpha

smis_jacobian
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o smis_jacobian_k2_stratified.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s smis_jacobian -k 2 -x

smis_all
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o smis_all_k2_stratified.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s smis_all -k 2 -x

cmis
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o cmis.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s cmis

average
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o Average.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s average

proxy
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o proxy.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s proxy_sample

proxy-sum
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o proxy_sum.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s proxy_sample -x

cargo run --features="pbrt progress-bar" --release --example=cli -- -t -1 -n 1 -o proxy.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 1 -s proxy_sample >a.txt

sims Elapsed Integrator: 27328 ms 2
smis Elapsed Integrator: 31123 ms 3



