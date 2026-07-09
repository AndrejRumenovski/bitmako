# Deploying a public BitMako demo

Step-by-step for building a 20M-compound subset of the full 1.36B corpus and
hosting it as a public, clickable demo at **bitmako.duckdns.org**, at $0 cost.
Steps 1 runs on your workstation (where the real `data/` files live); steps
2+ run on a free Oracle Cloud "Always Free" VM.

Stack: **Oracle Cloud Free Tier** (Ampere A1, permanently free — 4 OCPU/24GB
RAM/200GB disk, comfortably fits the ~4.5GB subset) + **DuckDNS** (free
dynamic DNS, gives you `bitmako.duckdns.org` pointing at the VM) + **nginx +
certbot** for real TLS on that free subdomain.

## 1. Build the 20M-compound subset (on your workstation)

Resulting files: ~1 GB index, ~2.5 GB fp-store, ~400 MB prop-store, ~500 MB
skip index — about 4.5 GB total.

```powershell
cargo build --release

# Carve the subset out of the real compounds.lance (no re-ingest needed).
.\target\release\bitmako.exe extract-subset `
  --lance data\compounds.lance --output data\demo\compounds.lance --limit 20000000

# Build the same four artifacts you'd build for the full corpus, just against
# the subset — identical commands, smaller/faster.
.\target\release\bitmako.exe build-index `
  --lance data\demo\compounds.lance --output data\demo\compounds.bitmako
.\target\release\bitmako.exe build-skip `
  --index data\demo\compounds.bitmako --output data\demo\compounds.skip
.\target\release\bitmako.exe build-fp-store `
  --lance data\demo\compounds.lance --output data\demo\compounds.fp
.\target\release\bitmako.exe build-prop-store `
  --lance data\demo\compounds.lance --output data\demo\compounds.prop
```

**Test it locally first** before shipping anything to a server:

```powershell
.\target\release\bitmako.exe serve `
  --index data\demo\compounds.bitmako --skip data\demo\compounds.skip `
  --fp-store data\demo\compounds.fp --prop-store data\demo\compounds.prop `
  --lance data\demo\compounds.lance --port 8080 `
  --demo-notice "Demo instance — a subset of the full 1.36B-compound corpus."
```

Visit `http://localhost:8080/` — confirm search works and the amber demo
banner shows up.

**Known characteristic of a `--limit`-only (no `--stride`) subset:** it's the
*first* N rows of the catalog, not a spread-out sample — confirmed against
the actual 20M build that this means some generic queries (plain aspirin,
benzene) return few or zero hits at realistic thresholds, since this slice
happens to be a narrow batch of the catalog rather than the full chemical
diversity. This isn't a bug (verified via a self-match search: a SMILES known
to be in the subset returns a correct `score=1.0` exact match, and the WAND
vs. brute-force correctness tests in `tests/search_correctness.rs` pass
independent of which compounds are in the corpus). For a demo query that
reliably returns a rich, populated result set in *this* 20M build, try:

```
CC(=O)OCCC(=O)NC1=CC(C)=NO1
```

(10 hits from score 1.0 down to 0.6 at threshold 0.5 — already reflected in
the `--demo-notice` text in `bitmako-demo.service`.) If you rebuild the subset
with different `--limit`/`--stride` values, re-derive a good example query the
same way: run a few searches locally and pick one with a populated neighborhood
before shipping.

## 2. Create the free Oracle Cloud VM

1. Sign up at [oracle.com/cloud/free](https://www.oracle.com/cloud/free/) (needs a card for identity verification, but the Always Free tier is never billed unless you explicitly upgrade — no auto-charge after a trial).
2. **Compute → Instances → Create instance.**
   - Image: **Ubuntu 22.04** (or the latest LTS offered).
   - Shape: click "Change shape" → **Ampere (VM.Standard.A1.Flex)** → this is the *Always Free* shape. Set 4 OCPUs / 24 GB RAM (the free-tier maximum).
   - Boot volume: bump to 200 GB (still free-tier eligible).
   - Add your SSH public key (generate one with `ssh-keygen` if you don't have one).
3. **Note the instance's public IP** once it's running.
4. **Open the firewall** — Oracle blocks ports by default at the network level, separate from the OS firewall:
   - Go to the instance's subnet → **Security Lists** → default security list → **Add Ingress Rules**: allow TCP 80 and TCP 443 from `0.0.0.0/0`.
   - On the VM itself, Ubuntu's `iptables`/`netfilter-persistent` also blocks by default on Oracle images — the commands in step 3 below open it there too.

## 3. Provision the VM

```bash
ssh ubuntu@YOUR_VM_IP

sudo apt update && sudo apt install -y build-essential pkg-config libssl-dev \
  protobuf-compiler git nginx certbot python3-certbot-nginx

# Oracle's Ubuntu images firewall everything by default at the OS level too.
sudo iptables -I INPUT -p tcp --dport 80 -j ACCEPT
sudo iptables -I INPUT -p tcp --dport 443 -j ACCEPT
sudo netfilter-persistent save 2>/dev/null || true

curl https://sh.rustup.rs -sSf | sh -s -- -y
source "$HOME/.cargo/env"

sudo useradd -m -s /usr/sbin/nologin bitmako
sudo mkdir -p /opt/bitmako/{bin,data/demo}
sudo chown -R bitmako:bitmako /opt/bitmako

git clone https://github.com/AndrejRumenovski/bitmako.git ~/bitmako-src
cd ~/bitmako-src
cargo build --release
sudo cp target/release/bitmako /opt/bitmako/bin/
```

## 4. Copy the subset data up

From your workstation:

```powershell
scp -r data\demo\* ubuntu@YOUR_VM_IP:/tmp/bitmako-demo-data/
```

Then on the VM:

```bash
sudo mv /tmp/bitmako-demo-data/* /opt/bitmako/data/demo/
sudo chown -R bitmako:bitmako /opt/bitmako/data
```

## 5. Install the systemd service

```bash
sudo cp ~/bitmako-src/deploy/bitmako-demo.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable --now bitmako-demo
sudo systemctl status bitmako-demo   # should show "active (running)"
curl http://127.0.0.1:8080/health    # sanity check before exposing publicly
```

## 6. Set up DuckDNS for bitmako.duckdns.org

1. Go to [duckdns.org](https://www.duckdns.org/), sign in (GitHub/Google login), and claim the subdomain **bitmako** — this gives you `bitmako.duckdns.org` pointing wherever you tell it, for free, forever.
2. On the DuckDNS dashboard, set the IP field to your Oracle VM's public IP (or use their update URL — see below to automate it).
3. Keep it pointed at the VM automatically (the VM's IP can change on reboot depending on your Oracle networking setup — a reserved/static public IP avoids this, which Oracle's free tier supports; assign one under **Networking → Reserved Public IPs**). If you don't reserve one, set up a cron job on the VM to keep DuckDNS updated:

```bash
mkdir -p ~/duckdns
cat > ~/duckdns/update.sh <<'EOF'
echo url="https://www.duckdns.org/update?domains=bitmako&token=YOUR_DUCKDNS_TOKEN&ip=" | curl -k -o ~/duckdns/duck.log -K -
EOF
chmod +x ~/duckdns/update.sh
(crontab -l 2>/dev/null; echo "*/5 * * * * ~/duckdns/update.sh >/dev/null 2>&1") | crontab -
```

(Get `YOUR_DUCKDNS_TOKEN` from the DuckDNS dashboard.)

## 7. Put nginx + TLS in front

`deploy/nginx-bitmako.conf` is already configured for `bitmako.duckdns.org`:

```bash
sudo cp ~/bitmako-src/deploy/nginx-bitmako.conf /etc/nginx/sites-available/bitmako-demo
sudo ln -s /etc/nginx/sites-available/bitmako-demo /etc/nginx/sites-enabled/
sudo nginx -t && sudo systemctl reload nginx

sudo certbot --nginx -d bitmako.duckdns.org   # provisions + auto-renews TLS
```

Visit **https://bitmako.duckdns.org/** — you should see the live search UI.

## Note on rate limiting in this topology

BitMako's own per-IP limiter (`src/api.rs`) identifies clients by TCP peer
address. Once nginx is the reverse proxy, every request's peer address is
nginx itself (`127.0.0.1`) — so BitMako's limiter becomes one **shared**
budget across all visitors combined (30 req/min total), not 30-per-visitor.
The nginx config (`limit_req_zone $binary_remote_addr`) is what actually
enforces per-visitor limits in this deployment; it sees the real client IP.
Both layers still help (nginx = per-visitor, BitMako = a global ceiling), just
know they're not doing the same job.

## Updating the demo later

```bash
cd ~/bitmako-src && git pull && cargo build --release
sudo cp target/release/bitmako /opt/bitmako/bin/
sudo systemctl restart bitmako-demo
```
